#' @name gl.report.replicates.KING
#' @title Identifies duplicate samples by KING-robust kinship
#' @family matched report
#'
#' @description
#' Scans a genlight object for duplicate samples (two libraries of one
#' animal) by KING-robust kinship and the heterozygote share, computed from
#' genotype likelihoods, read counts or called genotypes, and suggests which
#' sample of each group of duplicates to remove. Unlike the identity screen
#' of \code{\link{gl.report.replicates}}, it also finds duplicates in which
#' one sample is poor: missed heterozygotes in a poor sample lower the share
#' of identical genotypes, but much less so the kinship estimated from
#' likelihoods or the heterozygote share.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param input Data the kinship is computed from: "likelihoods" (the
#' genotype likelihoods kept by \code{gl.read.vcf(likelihoods = TRUE)}),
#' "counts" (the read counts kept by \code{gl.read.vcf(read.counts = TRUE)}),
#' "calls" (the called genotypes) or "auto", the first of these the object
#' holds [default "auto"].
#' @param error_rate Sequencing error rate of the binomial likelihoods built
#' from read counts; used only with read counts [default 0.001].
#' @param maf Minimum minor allele frequency of the loci used [default 0.05].
#' @param loc_threshold Minimum number of loci two samples must share for
#' their pair to be assessed [default 100].
#' @param king_threshold KING-robust kinship above which a pair is a
#' duplicate [default 0.354].
#' @param het_threshold Heterozygote share above which a pair is a duplicate
#' [default 0.80].
#' @param plot.out If TRUE, the plot is produced and returned; it is also
#' printed unless verbose is 0 [default TRUE].
#' @param plot_theme User specified theme [default theme_dartR()].
#' @param plot_colors Vector with two colour names: flagged pairs and other
#' pairs in the scatter plot, borders and fill of the histogram
#' [default c("#2171B5", "#6BAED6")].
#' @param bins Number of bins of the histogram [default 100].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' \strong{Pair statistics}
#'
#' At each locus the genotypes of two samples (0, 1 or 2 copies of the ALT
#' allele) fall into nine combinations, A to I (Waples et al. 2019; rows are
#' the first sample, columns the second):
#' \tabular{lccc}{
#'            \tab 0 \tab 1 \tab 2 \cr
#'   0 \tab A \tab D \tab G \cr
#'   1 \tab B \tab E \tab H \cr
#'   2 \tab C \tab F \tab I
#' }
#' From their counts over the loci the two samples share, the function
#' reports KING-robust kinship, (E - 2(C + G)) / (B + D + H + F + 2E)
#' (Manichaikul et al. 2010); R0 = (C + G) / E and R1 = E / (B + C + D + F +
#' G + H) (Waples et al. 2019); and the heterozygote share, E / min(B + E +
#' H, D + E + F): among the loci where the sample with fewer heterozygous
#' calls is heterozygous, the share where the other is heterozygous too.
#' None of them needs allele frequencies.
#'
#' With called genotypes the counts are tallied directly. With genotype
#' likelihoods or read counts they are the expected counts of the
#' two-dimensional site frequency spectrum of the pair, estimated by EM from
#' the likelihoods as in NgsRelate (Hangh\enc{ø}{o}j et al. 2019), so that a
#' heterozygote missed at low depth counts as uncertain rather than as a
#' homozygote. Likelihoods are 10^(-PL/10); from read counts R (REF) and A
#' (ALT) they are (1 - e)^R e^A, 0.5^(R + A) and e^R (1 - e)^A for 0/0, 0/1
#' and 1/1, with e = \code{error_rate}. The EM runs in C++ (Rcpp and
#' RcppParallel), compiled the first time the function is called in a
#' session, which takes a little while.
#'
#' Only loci with a minor allele frequency of at least \code{maf} are used:
#' from genotype likelihoods the frequency is estimated by EM, as in
#' NgsRelate; otherwise it is taken from the called genotypes. Pairs that
#' share fewer than \code{loc_threshold} of these loci get NA and are never
#' flagged.
#'
#' \strong{Duplicates}
#'
#' A pair is flagged as a duplicate when its KING-robust kinship exceeds
#' \code{king_threshold} or its heterozygote share exceeds
#' \code{het_threshold}. The heterozygote share guards against heterozygote
#' dropout, which turns heterozygotes into homozygotes and rarely the
#' reverse: conditioning on the heterozygous calls of the poorer sample
#' cancels its dropout. It is expected to be about 1 for duplicates, 0.5 for
#' parent and offspring, at most 0.625 for full siblings, and 2pq for
#' unrelated pairs. From called genotypes, KING-robust kinship underestimates
#' duplicates in which a sample has heterozygote dropout; the heterozygote
#' share corrects for dropout in the poorer sample of a pair, but not when
#' both samples are poor, so use genotype likelihoods where possible.
#'
#' Flagged pairs are linked into groups (connected components). In each
#' group the sample with the lowest proportion of missing genotypes is kept
#' (ties keep the name that sorts first) and the others are suggested for
#' removal. A group can also link two different animals through a sample
#' that is a mixture of them (contamination, see
#' \code{\link{gl.report.contamination}}); the function warns and lists the
#' pairs in a group that are not flagged themselves, which should be checked
#' before samples are removed.
#'
#' \strong{Kinship classes}
#'
#' Each pair gets a class from the KING bins of Manichaikul et al. (2010):
#' first degree above 0.177, second degree above 0.0884, third degree above
#' 0.0442 and unrelated otherwise; flagged pairs are "duplicate", including
#' those flagged only by the heterozygote share. KING-robust kinship assumes
#' a single population: between divergent lineages it is strongly negative
#' (values down to -3 are expected), which still means unrelated.
#'
#' \strong{Plot}
#'
#' The top panel plots every pair by KING-robust kinship and heterozygote
#' share, with dashed lines at \code{king_threshold}, at 0.177 and at
#' \code{het_threshold}; flagged pairs are labelled when there are 30 or
#' fewer and package ggrepel is installed. The kinship axis is cut at -1,
#' and a caption gives the number of pairs below it. Duplicates sit to the
#' right (from likelihoods) or at the top (from calls with dropout). The
#' bottom panel is a histogram of the statistic that separates duplicates:
#' the kinship of the pairs above 0, with the bounds of the classes, or,
#' from called genotypes, the heterozygote share of the pairs above 0.4,
#' with lines at 0.5 (parent and offspring), 0.625 (the largest expected
#' for full siblings) and \code{het_threshold}. A gap between the relatives
#' and the duplicates means the duplicate calls are clear-cut.
#'
#' @return A list, returned invisibly, with:
#' \itemize{
#' \item table.rep: the flagged pairs, with their statistics, class, group,
#' the proportion of missing genotypes of each sample (ind1_miss,
#' ind2_miss) and ind_to_drop, the sample of the pair with more missing
#' genotypes (ties: the name that sorts last). In a group where every pair
#' is flagged, ind_to_drop lists exactly ind.list.drop.
#' \item ind.list.drop: the samples suggested for removal.
#' \item ind.list.rep: for each sample in a flagged pair, its duplicates.
#' \item table.all: every pair (ind1 before ind2 in the object), with nloc
#' (shared loci), KING, R0, R1, het_share, class and flagged.
#' \item ind.miss: the proportion of missing genotypes of each sample.
#' \item input: the data used, "likelihoods", "counts" or "calls".
#' \item parameters: the thresholds used and n.sites, the number of loci
#' that passed \code{maf}.
#' \item plot: the plot (a patchwork of two ggplots), or NULL when
#' \code{plot.out = FALSE}.
#' }
#'
#' @references
#' \itemize{
#' \item Hangh\enc{ø}{o}j K, Moltke I, Andersen PA, Manica A, Korneliussen TS
#' (2019) Fast and accurate relatedness estimation from high-throughput
#' sequencing data in the presence of inbreeding. GigaScience 8: giz034.
#' \item Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM
#' (2010) Robust relationship inference in genome-wide association studies.
#' Bioinformatics 26: 2867-2873.
#' \item Waples RK, Albrechtsen A, Moltke I (2019) Allele frequency-free
#' inference of close familial relationships from genotypes or low-depth
#' sequencing data. Molecular Ecology 28: 35-48.
#' }
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' res <- gl.report.replicates.KING(platypus.gl)
#' head(res$table.all)
#' }
#'
#' @seealso \code{\link{gl.filter.replicates.KING}},
#' \code{\link{gl.report.replicates}}, \code{\link{gl.read.vcf}}
#' @export

gl.report.replicates.KING <- function(x,
                                      input = c("auto", "likelihoods",
                                                "counts", "calls"),
                                      error_rate = 0.001,
                                      maf = 0.05,
                                      loc_threshold = 100,
                                      king_threshold = 0.354,
                                      het_threshold = 0.80,
                                      plot.out = TRUE,
                                      plot_theme = theme_dartR(),
                                      plot_colors = c("#2171B5", "#6BAED6"),
                                      bins = 100,
                                      verbose = NULL) {

  # Pre-allocate variables to avoid "no visible binding" notes
  flagged <- group <- KING <- ind1_miss <- ind2_miss <- ind_to_drop <- NULL

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, build = "Jody", verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  input <- match.arg(input)
  number <- function(v) is.numeric(v) && length(v) == 1 && !is.na(v)
  if (!number(error_rate) || error_rate <= 0 || error_rate >= 0.5) {
    stop(error("  error_rate must lie in (0, 0.5)\n"))
  }
  if (!number(maf) || maf < 0 || maf >= 0.5) {
    stop(error("  maf must lie in [0, 0.5)\n"))
  }
  if (!number(loc_threshold) || loc_threshold < 1) {
    stop(error("  loc_threshold must be a number of at least 1\n"))
  }
  if (!number(king_threshold) || !number(het_threshold)) {
    stop(error("  king_threshold and het_threshold must be numbers\n"))
  }
  if (!number(bins) || bins < 1) {
    stop(error("  bins must be a positive number\n"))
  }
  if (length(plot_colors) < 2) {
    stop(error("  plot_colors must hold two colours\n"))
  }

  gm <- x@other$geno.metrics
  holds <- function(nm) {
    all(nm %in% names(gm)) &&
      all(vapply(gm[nm], function(m) identical(dim(m), c(nInd(x), nLoc(x))),
                 logical(1)))
  }
  pl.names <- c("pl.ref", "pl.het", "pl.alt")
  if (input == "auto") {
    input <- if (holds(pl.names)) {
      "likelihoods"
    } else if (holds(c("ref.count", "alt.count"))) {
      "counts"
    } else {
      "calls"
    }
  }
  if (input == "likelihoods" && !holds(pl.names)) {
    stop(error("  input = \"likelihoods\" needs the genotype likelihoods kept by gl.read.vcf(likelihoods = TRUE) in other$geno.metrics\n"))
  }
  if (input == "counts" && !holds(c("ref.count", "alt.count"))) {
    stop(error("  input = \"counts\" needs the read counts kept by gl.read.vcf(read.counts = TRUE) in other$geno.metrics\n"))
  }
  if (input != "calls") {
    for (pkg in c("Rcpp", "RcppParallel")) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(error("  Package", pkg, "needed for this function to work. Please install it.\n"))
      }
    }
  }

  if (verbose >= 2) {
    cat(report("  Kinship from", switch(input,
      likelihoods = "genotype likelihoods (other$geno.metrics: pl.ref, pl.het, pl.alt)",
      counts = paste0("read counts (other$geno.metrics: ref.count, alt.count; error rate ", error_rate, ")"),
      calls = "called genotypes"), "\n"))
  }
  if (input == "calls" && verbose >= 1) {
    cat(warn("  Warning: from called genotypes, KING-robust kinship underestimates duplicates in which a sample has missed heterozygotes (heterozygote dropout). The heterozygote share corrects for dropout in the poorer sample of a pair, but not when both samples are poor. Genotype likelihoods (gl.read.vcf(likelihoods = TRUE)) avoid this.\n"))
  }

  # DO THE JOB
  G <- as.matrix(x)
  ind.miss <- rowMeans(is.na(G))
  names(ind.miss) <- indNames(x)

  # Loci by minor allele frequency
  if (input == "likelihoods") {
    L <- lapply(gm[pl.names], function(m) t(10^(-m / 10)))
    f <- .king_freq(L)
  } else {
    f <- colMeans(G, na.rm = TRUE) / 2
  }
  use <- which(f >= maf & f <= 1 - maf)
  if (verbose >= 2) {
    cat(report("  Using", length(use), "of", nLoc(x),
               "loci with a minor allele frequency of at least", maf, "\n"))
  }

  # One row per unordered pair, ind1 earlier in the object than ind2
  n <- nInd(x)
  pairs <- which(upper.tri(matrix(TRUE, n, n)), arr.ind = TRUE)
  if (input == "calls") {
    s9 <- .king_calls9(G, use, pairs)
    nloc <- rowSums(s9)
  } else {
    if (input == "likelihoods") {
      L <- lapply(L, function(m) m[use, , drop = FALSE])
    } else {
      L <- .king_counts_lik(gm, use, error_rate)
    }
    has <- !is.na(L[[1]]) & !is.na(L[[2]]) & !is.na(L[[3]])
    L <- lapply(L, function(m) {
      m[!has] <- NA
      m
    })
    nloc <- crossprod(has * 1)[pairs]
    s9 <- matrix(NA_real_, nrow(pairs), 9, dimnames = list(NULL, LETTERS[1:9]))
    run <- which(nloc >= loc_threshold)
    if (length(run) > 0) {
      if (verbose >= 2) {
        cat(report("  Estimating the genotype combinations of", length(run),
                   "pairs by EM\n"))
      }
      em <- .king_em(L, pairs[run, , drop = FALSE] - 1L)
      s9[run, ] <- em[, 2:10]
      if (any(em[, 11] > 5000) && verbose >= 1) {
        cat(warn("  Warning: the EM did not converge in 5000 iterations for",
                 sum(em[, 11] > 5000), "pairs\n"))
      }
    }
  }
  st <- .king_stats(s9)
  st[nloc < loc_threshold, ] <- NA

  tab <- data.table(ind1 = indNames(x)[pairs[, 1]],
                    ind2 = indNames(x)[pairs[, 2]],
                    nloc = as.numeric(nloc))
  tab <- cbind(tab, st)
  tab[, flagged := .king_flag(tab, king_threshold, het_threshold,
                              loc_threshold)]
  tab[, class := .king_class(KING, flagged)]
  setcolorder(tab, c("ind1", "ind2", "nloc", "KING", "R0", "R1", "het_share",
                     "class", "flagged"))

  # Groups of duplicates and the samples to drop
  grp <- .king_groups(tab$ind1, tab$ind2, tab$flagged, ind.miss)
  rep.tab <- tab[tab$flagged, !"flagged"]
  rep.tab[, group := grp$group[tab$flagged]]
  rep.tab[, ind1_miss := unname(ind.miss[rep.tab$ind1])]
  rep.tab[, ind2_miss := unname(ind.miss[rep.tab$ind2])]
  rank1 <- unname(grp$rank[rep.tab$ind1])
  rank2 <- unname(grp$rank[rep.tab$ind2])
  rep.tab[, ind_to_drop := ifelse(rank1 > rank2, rep.tab$ind1, rep.tab$ind2)]
  rep.tab <- rep.tab[order(group, -KING)]

  members <- indNames(x)[indNames(x) %in% c(rep.tab$ind1, rep.tab$ind2)]
  ind.list.rep <- lapply(members, function(id) {
    partner <- c(rep.tab$ind2[rep.tab$ind1 == id], rep.tab$ind1[rep.tab$ind2 == id])
    indNames(x)[indNames(x) %in% partner]
  })
  names(ind.list.rep) <- members

  .king_chain_warn(tab, grp, verbose)

  # PLOT
  p <- NULL
  if (plot.out) {
    p <- .king_plot(tab, input, king_threshold, het_threshold, plot_theme,
                    plot_colors, bins)
    if (verbose >= 1) {
      suppressWarnings(print(p))
    }
  }

  # REPORT THE RESULTS
  if (verbose >= 3) {
    assessed <- !is.na(tab$KING)
    by.king <- (tab$KING > king_threshold) %in% TRUE
    cat("  Pairs compared:", nrow(tab), "\n")
    cat("  Pairs sharing at least", loc_threshold, "loci:", sum(assessed), "\n")
    cat("  Pairs flagged as duplicates:", nrow(rep.tab), "(",
        sum(by.king), "by KING,", sum(tab$flagged & !by.king),
        "by the heterozygote share only )\n")
    cat("  Groups of duplicates:", length(unique(rep.tab$group)), "\n")
    cat(report("  Samples suggested for removal:", length(grp$drop), "\n"))
    if (nrow(rep.tab) > 0) {
      print(rep.tab)
    }
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  invisible(list(
    table.rep = rep.tab[],
    ind.list.drop = grp$drop,
    ind.list.rep = ind.list.rep,
    table.all = tab[],
    ind.miss = ind.miss,
    input = input,
    parameters = list(error_rate = if (input == "counts") error_rate else NA,
                      maf = maf,
                      loc_threshold = loc_threshold,
                      king_threshold = king_threshold,
                      het_threshold = het_threshold,
                      n.sites = length(use)),
    plot = p
  ))
}

## Internal helpers of gl.report.replicates.KING and
## gl.filter.replicates.KING

## Allele frequency of each site by EM from genotype likelihoods with
## Hardy-Weinberg genotype priors, as in NgsRelate. L: the likelihoods of
## 0/0, 0/1 and 1/1 as three sites x individuals matrices, NA where an
## individual has no data.
.king_freq <- function(L, tol = 1e-8, maxit = 200) {
  f <- rep(0.25, nrow(L[[1]]))
  for (it in seq_len(maxit)) {
    w0 <- L[[1]] * (1 - f)^2
    w1 <- L[[2]] * 2 * f * (1 - f)
    w2 <- L[[3]] * f^2
    fn <- rowMeans((w1 + 2 * w2) / (w0 + w1 + w2), na.rm = TRUE) / 2
    done <- max(abs(fn - f), na.rm = TRUE) < tol
    f <- fn
    if (done) break
  }
  f
}

## Binomial genotype likelihoods from read counts (sites x individuals),
## computed on the log scale and scaled to a maximum of 1 so that deep
## genotypes do not underflow; NA where a genotype has no reads
.king_counts_lik <- function(gm, use, e) {
  R <- t(gm$ref.count[, use, drop = FALSE])
  A <- t(gm$alt.count[, use, drop = FALSE])
  l <- list(R * log1p(-e) + A * log(e),
            (R + A) * log(0.5),
            R * log(e) + A * log1p(-e))
  top <- do.call(pmax, l)
  miss <- is.na(top) | (R + A) == 0
  lapply(l, function(m) {
    m <- exp(m - top)
    m[miss] <- NA
    m
  })
}

## Counts of the nine genotype combinations A-I of each pair from called
## genotypes, by matrix products over blocks of loci. D, G and H are the
## transposes of B, C and F.
.king_calls9 <- function(G, use, pairs, block = 5000L) {
  n <- nrow(G)
  acc <- rep(list(matrix(0, n, n)), 6)
  for (b in split(use, ceiling(seq_along(use) / block))) {
    g <- G[, b, drop = FALSE]
    ind <- lapply(0:2, function(v) {
      m <- (g == v) * 1
      m[is.na(m)] <- 0
      m
    })
    acc[[1]] <- acc[[1]] + tcrossprod(ind[[1]], ind[[1]])
    acc[[2]] <- acc[[2]] + tcrossprod(ind[[2]], ind[[1]])
    acc[[3]] <- acc[[3]] + tcrossprod(ind[[3]], ind[[1]])
    acc[[4]] <- acc[[4]] + tcrossprod(ind[[2]], ind[[2]])
    acc[[5]] <- acc[[5]] + tcrossprod(ind[[3]], ind[[2]])
    acc[[6]] <- acc[[6]] + tcrossprod(ind[[3]], ind[[3]])
  }
  ji <- pairs[, 2:1, drop = FALSE]
  cbind(A = acc[[1]][pairs], B = acc[[2]][pairs], C = acc[[3]][pairs],
        D = acc[[2]][ji], E = acc[[4]][pairs], F = acc[[5]][pairs],
        G = acc[[3]][ji], H = acc[[5]][ji], I = acc[[6]][pairs])
}

## Expected counts of the nine genotype combinations of each pair by EM from
## genotype likelihoods (Waples et al. 2019; NgsRelate), in C++. pairs:
## 0-based individuals. Returns pairs x 11: shared sites, A-I, iterations
## (above maxit when the EM did not converge).
.king_em <- function(L, pairs, tol = 1e-9, maxit = 5000L) {
  # Stub to satisfy package checks when defining the real function below
  kingEMParallel <- function() {}
  suppressWarnings(
    Rcpp::cppFunction(
      code = '
      #include <Rcpp.h>
      #include <RcppParallel.h>
      #include <algorithm>
      #include <cmath>
      #include <vector>
      using namespace Rcpp;
      using namespace RcppParallel;

      struct KingEM : public Worker {
        const RMatrix<double> L0, L1, L2;  // sites x individuals
        const RVector<int> p1, p2;          // 0-based individuals of each pair
        const double tol;
        const int maxit;
        RMatrix<double> out;

        KingEM(const NumericMatrix& L0_, const NumericMatrix& L1_,
               const NumericMatrix& L2_, const IntegerVector& p1_,
               const IntegerVector& p2_, double tol_, int maxit_,
               NumericMatrix& out_)
          : L0(L0_), L1(L1_), L2(L2_), p1(p1_), p2(p2_), tol(tol_),
            maxit(maxit_), out(out_) {}

        void operator()(std::size_t begin, std::size_t end) override {
          const std::size_t ns = L0.nrow();
          std::vector<double> M;
          M.reserve(ns * 9);
          for (std::size_t q = begin; q < end; ++q) {
            const std::size_t i = p1[q], j = p2[q];
            // likelihood of each combination at each shared site, in the
            // order A-I: individual 1 varies fastest
            M.clear();
            for (std::size_t s = 0; s < ns; ++s) {
              if (std::isnan(L0(s, i)) || std::isnan(L0(s, j))) continue;
              const double a[3] = {L0(s, i), L1(s, i), L2(s, i)};
              const double b[3] = {L0(s, j), L1(s, j), L2(s, j)};
              for (int g2 = 0; g2 < 3; ++g2)
                for (int g1 = 0; g1 < 3; ++g1)
                  M.push_back(a[g1] * b[g2]);
            }
            const std::size_t n = M.size() / 9;
            double p[9], acc[9], pn[9];
            for (int k = 0; k < 9; ++k) p[k] = 1.0 / 9.0;
            int it = 0;
            if (n > 0) {
              // p is replaced by the mean posterior of each combination
              for (it = 1; it <= maxit; ++it) {
                for (int k = 0; k < 9; ++k) acc[k] = 0.0;
                for (std::size_t s = 0; s < n; ++s) {
                  const double* m = &M[s * 9];
                  double d = 0.0;
                  for (int k = 0; k < 9; ++k) d += m[k] * p[k];
                  for (int k = 0; k < 9; ++k) acc[k] += m[k] / d;
                }
                double diff = 0.0;
                for (int k = 0; k < 9; ++k) {
                  pn[k] = p[k] * (acc[k] / n);
                  diff = std::max(diff, std::fabs(pn[k] - p[k]));
                }
                for (int k = 0; k < 9; ++k) p[k] = pn[k];
                if (diff < tol) break;
              }
            }
            out(q, 0) = n;
            for (int k = 0; k < 9; ++k) out(q, k + 1) = p[k] * n;
            out(q, 10) = it;
          }
        }
      };

      // [[Rcpp::export]]
      NumericMatrix kingEMParallel(const NumericMatrix& L0,
                                   const NumericMatrix& L1,
                                   const NumericMatrix& L2,
                                   const IntegerVector& p1,
                                   const IntegerVector& p2,
                                   double tol, int maxit) {
        NumericMatrix out(p1.size(), 11);
        KingEM worker(L0, L1, L2, p1, p2, tol, maxit, out);
        parallelFor(0, p1.size(), worker);
        return out;
      }
      ',
      depends = "RcppParallel",
      plugin = "cpp11"
    )
  )
  kingEMParallel(L[[1]], L[[2]], L[[3]], as.integer(pairs[, 1]),
                 as.integer(pairs[, 2]), tol, as.integer(maxit))
}

## KING-robust kinship, R0, R1 and the heterozygote share from the nine
## counts (a matrix with columns A-I)
.king_stats <- function(s) {
  A <- s[, "A"]; B <- s[, "B"]; C <- s[, "C"]; D <- s[, "D"]; E <- s[, "E"]
  F <- s[, "F"]; G <- s[, "G"]; H <- s[, "H"]
  data.frame(KING = (E - 2 * (C + G)) / (B + D + H + F + 2 * E),
             R0 = (C + G) / E,
             R1 = E / (B + C + D + F + G + H),
             het_share = E / pmin(B + E + H, D + E + F))
}

## Pairs flagged as duplicates: enough shared loci and KING or the
## heterozygote share above its threshold
.king_flag <- function(tab, king_threshold, het_threshold, loc_threshold) {
  (tab$nloc >= loc_threshold &
     (tab$KING > king_threshold | tab$het_share > het_threshold)) %in% TRUE
}

## KING class of each pair (Manichaikul et al. 2010); flagged pairs are
## duplicates
.king_class <- function(king, flagged) {
  cl <- as.character(cut(king, c(-Inf, 0.0442, 0.0884, 0.177, Inf),
                         labels = c("unrelated", "third degree",
                                    "second degree", "first degree")))
  cl[flagged] <- "duplicate"
  cl
}

## Links flagged pairs into groups (connected components), keeps in each
## group the sample with the fewest missing genotypes (ties: the name that
## sorts first) and drops the others. Returns the group of every pair (NA
## outside groups), the samples kept and dropped, the rank of each sample in
## that order and the pairs inside a group that are not flagged.
.king_groups <- function(ind1, ind2, flagged, ind.miss) {
  ids <- names(ind.miss)
  rank <- integer(length(ids))
  rank[order(unname(ind.miss), ids, method = "radix")] <- seq_along(ids)
  names(rank) <- ids
  res <- list(group = rep(NA_integer_, length(ind1)), keep = character(0),
              drop = character(0), rank = rank, chain = integer(0))
  a <- ind1[flagged]
  b <- ind2[flagged]
  if (length(a) == 0) {
    return(res)
  }
  # union-find over the samples in flagged pairs, in object order, so that
  # each group is named after its earliest sample
  members <- ids[ids %in% c(a, b)]
  parent <- seq_along(members)
  find <- function(k) {
    while (parent[k] != k) k <- parent[k]
    k
  }
  ia <- match(a, members)
  ib <- match(b, members)
  for (q in seq_along(ia)) {
    ra <- find(ia[q])
    rb <- find(ib[q])
    if (ra != rb) parent[max(ra, rb)] <- min(ra, rb)
  }
  root <- vapply(seq_along(members), find, integer(1))
  gid <- stats::setNames(match(root, unique(root)), members)
  keep <- vapply(split(members, gid[members]),
                 function(m) m[which.min(rank[m])], character(1))
  g1 <- gid[ind1]
  g2 <- gid[ind2]
  same <- !is.na(g1) & !is.na(g2) & g1 == g2
  res$group <- unname(ifelse(same, g1, NA_integer_))
  res$keep <- unname(keep)
  res$drop <- setdiff(members, keep)
  res$chain <- which(same & !flagged)
  res
}

## Warns about pairs inside a group of duplicates that are not flagged
## themselves
.king_chain_warn <- function(tab, grp, verbose) {
  if (length(grp$chain) > 0 && verbose >= 1) {
    cat(warn("  Warning: these pairs are in the same group of duplicates but are not flagged themselves; they are linked only through other samples. A sample that is a mixture of two animals (contamination) can link them, so check them before removing samples:\n"))
    show <- data.frame(group = grp$group[grp$chain], ind1 = tab$ind1[grp$chain],
                       ind2 = tab$ind2[grp$chain], KING = tab$KING[grp$chain],
                       het_share = tab$het_share[grp$chain])
    print(show, row.names = FALSE, digits = 3)
  }
}

## Scatter of KING against the heterozygote share over a histogram of the
## statistic that separates duplicates: KING above 0, or from called
## genotypes, where dropout lowers KING, the heterozygote share above 0.4
.king_plot <- function(tab, input, king_threshold, het_threshold, plot_theme,
                       plot_colors, bins) {
  KING <- het_share <- status <- label <- value <- NULL
  d <- as.data.frame(tab[!is.na(tab$KING) & !is.na(tab$het_share), ])
  d$status <- factor(ifelse(d$flagged, "duplicate", "other"),
                     levels = c("duplicate", "other"))
  d <- d[order(d$flagged), ]
  # pairs between divergent lineages reach strongly negative values; the
  # axis is cut at -1 so that the duplicates and close relatives stay legible
  low <- max(min(c(d$KING, 0)), -1)
  hidden <- sum(d$KING < low)
  p1 <- ggplot(d, aes(KING, het_share, colour = status, size = status)) +
    geom_point(alpha = 0.7) +
    geom_vline(xintercept = c(0.177, king_threshold), linetype = "dashed") +
    geom_hline(yintercept = het_threshold, linetype = "dashed") +
    scale_colour_manual(values = c(duplicate = plot_colors[1],
                                   other = plot_colors[2]), drop = FALSE) +
    scale_size_manual(values = c(duplicate = 2.5, other = 1), drop = FALSE) +
    coord_cartesian(xlim = c(low, max(c(d$KING, king_threshold)))) +
    labs(x = "KING-robust kinship", y = "Heterozygote share", colour = NULL,
         size = NULL,
         caption = if (hidden > 0) {
           paste(hidden, "pairs with KING below -1 are not shown")
         }) +
    plot_theme
  fl <- d[d$flagged, ]
  if (nrow(fl) > 0 && nrow(fl) <= 30 &&
      requireNamespace("ggrepel", quietly = TRUE)) {
    fl$label <- paste(fl$ind1, fl$ind2, sep = " / ")
    # flagged pairs sit in the top right corner, so their labels go to the
    # empty area to their left, on a white background that keeps them
    # legible over the dashed lines
    p1 <- p1 + ggrepel::geom_label_repel(
      data = fl, aes(label = label), size = 3, show.legend = FALSE,
      fill = "#FFFFFFD9", label.size = NA, box.padding = 0.3,
      min.segment.length = 0, nudge_x = -0.3 * (max(d$KING) - low),
      ylim = c(min(fl$het_share, het_threshold) - 0.2, NA), force = 2,
      max.time = 2, max.overlaps = Inf, seed = 1)
  }
  if (input == "calls") {
    h <- data.frame(value = d$het_share[d$het_share > 0.4])
    cuts <- c(0.5, 0.625, het_threshold)
    band <- c("parent-offspring", "full siblings (max)", "duplicate")
    xlab <- "Heterozygote share (pairs above 0.4)"
  } else {
    h <- data.frame(value = d$KING[d$KING > 0])
    cuts <- c(0.0884, 0.177, king_threshold)
    band <- c("2nd degree", "1st degree", "duplicate")
    xlab <- "KING-robust kinship (pairs above 0)"
  }
  p2 <- ggplot(h, aes(value)) +
    (if (nrow(h) > 0) {
      geom_histogram(bins = bins, colour = plot_colors[1],
                     fill = plot_colors[2])
    }) +
    geom_vline(xintercept = cuts, linetype = "dashed") +
    annotate("text", x = cuts, y = Inf, vjust = 1.5, hjust = -0.05, size = 3,
             label = band) +
    labs(x = xlab, y = "Count") +
    plot_theme
  p1 / p2
}
