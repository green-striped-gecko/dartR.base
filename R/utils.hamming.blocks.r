#' @name utils.hamming.engine
#' @title Compiled Hamming comparison engine shared by gl.filter.hamming and
#' gl.report.hamming
#'
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' @details
#' Compiles (on first call in a session) and returns two Rcpp functions:
#' \itemize{
#' \item \code{dedup}: the block-hashing near-duplicate detector used by
#' \code{gl.filter.hamming}. Given a list of equal-length raw vectors ordered
#' from worst to best call rate and a mismatch threshold \code{k}, it returns
#' a logical keep vector where \code{FALSE} marks a locus with a kept partner
#' within \code{k} mismatches, plus a flag indicating whether the candidate
#' cap was hit.
#' \item \code{pairwise}: exact Hamming mismatch counts for a given matrix of
#' index pairs, used by \code{gl.report.hamming} for the distance
#' distribution.
#' }
#' The compiled functions are cached in a package-local environment, so the
#' C++ toolchain runs at most once per session. Requires the package
#' \code{Rcpp} and a working compiler toolchain; callers must guard with
#' \code{requireNamespace("Rcpp")}.
#'
#' @return A list with elements \code{dedup} and \code{pairwise}.
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal

.hamming_engine_cache <- new.env(parent = emptyenv())

utils.hamming.engine <- function() {
  if (!is.null(.hamming_engine_cache$engine)) {
    return(.hamming_engine_cache$engine)
  }

  dedup <- Rcpp::cppFunction(includes = '
#include <unordered_map>
#include <vector>
#include <cstdint>

// 64-bit FNV-1a hash over a raw slice [start, end)
static inline uint64_t fnv1a64(const uint8_t* p, int start, int end) {
  const uint64_t FNV_OFFSET = 1469598103934665603ULL;
  const uint64_t FNV_PRIME  = 1099511628211ULL;
  uint64_t h = FNV_OFFSET;
  for (int i = start; i < end; ++i) {
    h ^= (uint64_t)p[i];
    h *= FNV_PRIME;
  }
  return h;
}

// Exact Hamming mismatches with early stop at k (over [0, L))
static inline bool within_k_mism(const uint8_t* a, const uint8_t* b, int L, int k) {
  int mism = 0;
  for (int i = 0; i < L; ++i) {
    mism += (a[i] != b[i]);
    if (mism > k) return false;
  }
  return true;
}
', code = '
List filter_hamming_blocks_cpp(List raws_trimmed,
                              int k,
                              int max_candidates_cap = 5000) {
  const int n = raws_trimmed.size();
  LogicalVector keep(n, true);
  bool capped = false;

  if (n <= 1) return List::create(_["keep"] = keep, _["capped"] = capped);

  // Assume all sequences same length after trimming
  RawVector r0 = raws_trimmed[0];
  const int L = r0.size();
  if (L <= 0) return List::create(_["keep"] = keep, _["capped"] = capped);

  for (int i = 1; i < n; ++i) {
    RawVector ri = raws_trimmed[i];
    if (ri.size() != L) stop("All trimmed sequences must have identical length for this function.");
  }

  if (k < 0) k = 0;
  if (k >= L) {
    // Everything matches everything => keep only the last (best call rate)
    for (int i = 0; i < n - 1; ++i) keep[i] = false;
    return List::create(_["keep"] = keep, _["capped"] = capped);
  }

  const int B = k + 1;             // number of blocks
  std::vector<int> bstart(B), bend(B);
  {
    const int base = L / B;
    const int rem  = L % B;
    int s = 0;
    for (int b = 0; b < B; ++b) {
      const int len = base + (b < rem ? 1 : 0);
      bstart[b] = s;
      bend[b]   = s + len;
      s += len;
    }
  }

  // One hash table per block: hash -> vector of indices already inserted (later indices)
  std::vector< std::unordered_map<uint64_t, std::vector<int>> > tables(B);
  for (int b = 0; b < B; ++b) tables[b].reserve((size_t)n / 2);

  std::vector<int> cand;
  cand.reserve(1024);
  std::vector<int> seen(n, 0);
  int seen_token = 1;

  // Process from end to start: later indices are in the tables
  for (int i = n - 1; i >= 0; --i) {
    RawVector ai = raws_trimmed[i];
    const uint8_t* ap = (const uint8_t*)RAW(ai);

    cand.clear();
    ++seen_token;

    // Collect candidates from any matching block bucket
    for (int b = 0; b < B; ++b) {
      const uint64_t h = fnv1a64(ap, bstart[b], bend[b]);
      auto it = tables[b].find(h);
      if (it == tables[b].end()) continue;

      const std::vector<int>& bucket = it->second;
      for (int idx : bucket) {
        if (seen[idx] == seen_token) continue;
        seen[idx] = seen_token;
        cand.push_back(idx);
        if ((int)cand.size() >= max_candidates_cap) { capped = true; break; }
      }
      if ((int)cand.size() >= max_candidates_cap) break;
    }

    // Verify candidates with exact Hamming <= k (early exit in within_k_mism)
    bool found = false;
    for (int jj = 0; jj < (int)cand.size(); ++jj) {
      const int j = cand[jj];
      RawVector bj = raws_trimmed[j];
      const uint8_t* bp = (const uint8_t*)RAW(bj);
      if (within_k_mism(ap, bp, L, k)) {
        found = true;
        break;
      }
    }

    if (found) {
      keep[i] = false;
      continue; // do NOT insert deleted i
    }

    // Insert kept i into all block tables
    for (int b = 0; b < B; ++b) {
      const uint64_t h = fnv1a64(ap, bstart[b], bend[b]);
      tables[b][h].push_back(i);
    }
  }

  return List::create(_["keep"] = keep, _["capped"] = capped);
}
', depends = "Rcpp")

  pairwise <- Rcpp::cppFunction(code = '
IntegerVector pairwise_hamming_cpp(List raws_trimmed, IntegerMatrix pairs) {
  const int m = pairs.nrow();
  IntegerVector out(m);
  for (int p = 0; p < m; ++p) {
    RawVector a = raws_trimmed[pairs(p, 0) - 1];
    RawVector b = raws_trimmed[pairs(p, 1) - 1];
    const uint8_t* ap = (const uint8_t*)RAW(a);
    const uint8_t* bp = (const uint8_t*)RAW(b);
    const int L = a.size();
    int mism = 0;
    for (int i = 0; i < L; ++i) mism += (ap[i] != bp[i]);
    out[p] = mism;
  }
  return out;
}
', depends = "Rcpp")

  .hamming_engine_cache$engine <- list(dedup = dedup, pairwise = pairwise)
  .hamming_engine_cache$engine
}
