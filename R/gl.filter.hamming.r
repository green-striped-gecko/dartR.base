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

  # Compiled comparison engine shared with gl.report.hamming
  # (utils.hamming.engine compiles once per session and caches)
  filter_hamming_blocks_cpp <- utils.hamming.engine()$dedup


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
