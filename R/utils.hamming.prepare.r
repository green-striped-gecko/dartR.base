#' @name utils.hamming.prepare
#' @title Validates arguments and prepares trimmed sequences for
#' gl.filter.hamming and gl.report.hamming
#'
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' @details
#' Holds the one copy of the argument checks and sequence preparation shared
#' by \code{gl.filter.hamming} and \code{gl.report.hamming}. The report
#' simulates the filter exactly only if both prepare the sequences in the
#' same way, so neither function repeats this code.
#'
#' Each \code{TrimmedSequence} is upper-cased, the first \code{rs} bases are
#' skipped and the next \code{min.length} bases are kept. Loci whose
#' substring is exactly \code{min.length} bases long are comparable; they are
#' ordered from most to least missing data, so that the compiled engine,
#' which keeps the later of two duplicates, always keeps the locus with the
#' better call rate.
#'
#' @param x Name of the genlight object [required].
#' @param threshold Maximum number of mismatching bases [required].
#' @param rs Number of bases to skip from the start of the TrimmedSequence
#' [required].
#' @param min.length Length of the compared substring [required].
#'
#' @return A list with elements \code{raws} (raw vectors, one per locus),
#' \code{idx} (indices of comparable loci), \code{ord} (\code{idx} ordered
#' from worst to best call rate) and \code{n.short} (number of loci that are
#' not comparable).
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal

utils.hamming.prepare <- function(x, threshold, rs, min.length) {
  is.count <- function(v, min) {
    is.numeric(v) && length(v) == 1 && !is.na(v) && v >= min &&
      v == round(v)
  }

  if (length(x@other$loc.metrics$TrimmedSequence) == 0) {
    stop(error("Fatal Error: Data must include Trimmed Sequences\n"))
  }

  if (!is.count(rs, 0)) {
    stop(error(
      "Fatal Error: rs, the number of bases to skip (restriction site",
      "length), must be a whole number of 0 or more; usually it is less",
      "than 9\n"
    ))
  }

  if (!is.count(min.length, 1)) {
    stop(error(
      "Fatal Error: min.length, the number of bases compared, must be a",
      "whole number of 1 or more\n"
    ))
  }

  # Checked before the whole-number test so that legacy calls get the
  # migration message rather than a generic one
  if (is.numeric(threshold) && length(threshold) == 1 && !is.na(threshold) &&
      threshold > 0 && threshold < 1) {
    stop(error(
      "Fatal Error: threshold is the maximum number of mismatching bases",
      "(e.g. 3), not a proportion. Earlier versions of gl.filter.hamming",
      "took a proportion (default 0.2)\n"
    ))
  }

  if (!is.count(threshold, 0)) {
    stop(error(
      "Fatal Error: threshold must be a single whole number of mismatching",
      "bases, 0 or more\n"
    ))
  }

  # At threshold >= min.length every pair of comparable loci counts as a
  # duplicate, and all but one comparable locus would be removed
  if (threshold >= min.length) {
    stop(error(
      "Fatal Error: threshold (", threshold, ") must be smaller than",
      "min.length (", min.length, "); otherwise every comparable locus",
      "matches every other\n"
    ))
  }

  seqs <- toupper(as.character(x@other$loc.metrics$TrimmedSequence))
  trimmed <- substr(seqs, rs + 1, min.length + rs)
  raws <- lapply(trimmed, charToRaw)

  idx <- which(lengths(raws) == min.length)
  na.counts <- glNA(x)
  ord <- idx[order(na.counts[idx], decreasing = TRUE)]

  list(raws = raws,
       idx = idx,
       ord = ord,
       n.short = nLoc(x) - length(idx))
}
