#' @name gl.filter.hamming
#' @title Filters loci by trimmed-sequence similarity using Hamming distance
#' @family matched filter
#'
#' @description
#' Identifies loci with highly similar (near-duplicate) trimmed tag sequences and
#' removes redundant loci, preferentially retaining the locus with the better call
#' rate (fewer missing genotypes). Loci with near-identical tags are likely to be
#' fragments of the same genomic locus that the calling pipeline has split into
#' separate loci; retaining both yields pseudo-replicated, non-independent markers.
#' This function compares locus \code{TrimmedSequence} strings after skipping
#' a user-defined number of bases (the restriction site) and retaining a fixed-length
#' substring. Loci whose substrings are within \code{threshold} mismatches
#' (Hamming distance) are considered duplicates and one is dropped.
#'
#' @details
#' This function finds near-duplicate TrimmedSequences by first taking the same
#' fixed-length piece of sequence from every locus, then cutting that piece
#' into several smaller sections. It relies on a simple fact: if two sequences
#' differ by only a few letters (up to the threshold), then at least one of
#' those sections must be exactly the same in both. So it uses the sections as
#' 'signatures' to quickly shortlist only those loci that share an identical
#' section, instead of comparing every locus to every other one. For each
#' shortlisted pair, it then checks the two sequences letter-by-letter and
#' stops as soon as it can tell they differ by more than the allowed number.
#'
#' Before comparison, loci are ordered from most to least missing data (worst
#' to best call rate). When two loci are found to be within \code{threshold}
#' mismatches, the one with more missing genotypes is dropped, so the retained
#' locus of every duplicate pair is always the one with the better call rate.
#' Once a locus is dropped it is not used to match others.
#'
#' The function expects locus metrics to include \code{TrimmedSequence} in
#' \code{x@other$loc.metrics}.
#'
#' Only loci whose \code{TrimmedSequence} is long enough to yield a substring of
#' exactly \code{min.length} are compared. Loci with shorter sequences are not
#' compared and are retained.
#'
#' This function requires the package \code{Rcpp} and a working C++ compiler
#' toolchain (Rtools on Windows, Xcode Command Line Tools on macOS), as the
#' comparison engine is compiled at run time.
#'
#' @param x Name of the genlight object containing the SNP or SilicoDArT data
#'   [required].
#' @param threshold Maximum allowed Hamming distance (number of
#'   mismatching bases) between two trimmed sequences for them to be treated as
#'   duplicates; an integer >= 0, where 0 removes exact duplicates only.
#'   Note that in earlier versions of this function the threshold
#'   was a proportion of the sequence length (default 0.2); it is now a count
#'   of bases, consistent with \code{gl.report.hamming} [default 3].
#' @param rs Number of bases to skip from the start of the TrimmedSequence
#'  before extracting the comparison substring (i.e. restriction site length)
#'   [default 5].
#' @param min.length Length of the substring used for Hamming comparisons.
#'   Despite the name, this is an exact length, not a minimum: longer sequences
#'   are truncated to it, and only loci producing a substring of exactly this
#'   length are compared; others are retained [default 50].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  brief progress messages; 3, progress and results summary; 5, full report
#'   [default 2, unless specified using gl.set.verbosity].
#'
#' @return A \code{genlight} object with redundant loci removed.
#' @seealso \code{\link{gl.report.hamming}} to explore the distribution of
#'   sequence divergence before choosing a threshold.
#' @author Custodian: Arthur Georges -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # x must be a genlight with TrimmedSequence in x@other$loc.metrics
#' \donttest{
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' x <- platypus.gl
#' x2 <- gl.filter.hamming(x, threshold = 3, rs = 5, min.length = 50, verbose = 2)
#' }
#' @export

gl.filter.hamming <- function(x,
                              threshold = 3,
                              rs = 5,
                              min.length = 50,
                              verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   # build = "v.2023.3",
                   verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  if (!(requireNamespace("Rcpp", quietly = TRUE))) {
    stop(error(
      "Package Rcpp needed for this function to work. Please install it.\n"
    ))
  }

  if (length(x@other$loc.metrics$TrimmedSequence) == 0) {
    stop(error("Fatal Error: Data must include Trimmed Sequences\n"))
  }

  if (threshold < 0) {
    stop(error("Fatal Error: threshold must be a non-negative number of bases\n"))
  }

  if (threshold > 0 && threshold < 1) {
    stop(error(
      "Fatal Error: threshold is the maximum number of mismatching bases",
      " (e.g. 3), not a proportion. Earlier versions of this function took",
      " a proportion (default 0.2)\n"
    ))
  }

  # DO THE JOB

  n0 <- nLoc(x)

  #setup empty function for CRAN checks
  filter_hamming_blocks_cpp <- function() {  }

  Rcpp::cppFunction(includes = '
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


    seqs <- toupper(as.character(x@other$loc.metrics$TrimmedSequence))
    trimmed <- substr(seqs, rs + 1 , min.length + rs)
    raws <- lapply(trimmed, charToRaw)
    lens <- lengths(raws)

    idx <- which(lens == min.length)
    n.short <- n0 - length(idx)
    if (verbose >= 2 && n.short > 0) {
      cat(report(
        " ", n.short, "loci with a TrimmedSequence shorter than",
        min.length + rs, "bases were not compared and are retained\n"
      ))
    }

    # Order comparable loci from worst to best call rate. The C++ engine keeps
    # the later of two duplicates, so the better locus always survives and
    # remains available to match further duplicates.
    na.counts <- glNA(x)
    ord <- idx[order(na.counts[idx], decreasing = TRUE)]

    res <- filter_hamming_blocks_cpp(raws[ord], k = threshold,
                                     max_candidates_cap = 5000)

    if (res$capped && verbose >= 1) {
      cat(warn(
        "  Warning: more than 5000 candidate matches for at least one locus;",
        " some duplicates may have been missed\n"
      ))
    }

    drop.idx <- ord[!res$keep]

    if (length(drop.idx) > 0) {
      keep.idx <- setdiff(seq_len(n0), drop.idx)
      x2 <- x[, keep.idx]
      x2@other$loc.metrics <- x@other$loc.metrics[keep.idx, , drop = FALSE]
    } else {
      x2 <- x
    }

    # REPORT A SUMMARY
    if (verbose >= 3) {
      cat("\n  Summary of filtered dataset\n")
      cat(paste("    Initial No. of loci:", n0, "\n"))
      cat(paste("    Loci deleted", (n0 - nLoc(x2)), "\n"))
      cat(paste("    Final No. of loci:", nLoc(x2), "\n"))
      cat(paste("    No. of individuals:", nInd(x2), "\n"))
      cat(paste("    No. of populations: ", length(levels(factor(
        pop(x2)
      ))), "\n"))
    }

    # ADD TO HISTORY
    nh <- length(x2@other$history)
    x2@other$history[[nh + 1]] <- match.call()

    # FLAG SCRIPT END
    if (verbose > 0) {
      cat(report("Completed:", funname, "\n"))
    }

    # RETURN

    return(x2)
}
