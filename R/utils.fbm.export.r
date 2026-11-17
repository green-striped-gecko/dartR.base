## Exported block-wise helpers for file-backed (FBM) dartR objects, for use
## by the other dartRverse packages (which must not reach dartR.base
## internals with :::). Each is a thin wrapper of the internal helper in
## utils.fbm.r, so dartR.base and the downstream packages share one
## implementation.

#' @name utils.fbm.has
#' @title Tests whether a dartR object holds its genotypes in a file-backed
#' matrix
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' TRUE if \code{x} has an \code{fbm} slot holding a
#' \code{bigstatsr::FBM.code256}, FALSE for a gen-backed object (also for a
#' plain adegenet genlight, which has no such slot).
#' @param x A genlight or dartR object [required].
#' @return TRUE or FALSE.
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.fbm.has <- function(x) {
  .has_fbm(x)
}

#' @name utils.fbm.colwise
#' @title Applies a per-locus statistic in blocks of loci
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' Returns the same value as \code{FUN(as.matrix(x)[rows, ])}, where FUN
#' computes one value per locus (e.g. \code{colMeans(m, na.rm = TRUE)}). A
#' file-backed object is read in blocks of loci (size set by
#' \code{options(bigstatsr.block.sizeGB)}) and the results are concatenated,
#' so the whole genotype matrix is never decoded; blocks are decoded exactly
#' as by \code{as.matrix()} (integer storage, indNames/locNames as dimnames),
#' so the values are identical. A gen-backed object is decoded once. FUN may
#' return a matrix with one column per locus; the blocks are then bound with
#' \code{cbind()}.
#' @param x A dartR object [required].
#' @param FUN Function of a genotype matrix (individuals x loci) returning
#' one value (or one column) per locus [required].
#' @param rows Row indices (e.g. the individuals of one population), or NULL
#' for all individuals [default NULL].
#' @return A vector with one value per locus, or a matrix with one column
#' per locus.
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.fbm.colwise <- function(x, FUN, rows = NULL) {
  .gl_colwise(x, FUN, rows = rows)
}

#' @name utils.fbm.colblocks
#' @title Applies a function to blocks of loci and returns the per-block
#' results
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' As \code{\link{utils.fbm.colwise}}, but returns the list of per-block
#' results, for callers that combine them themselves (e.g. FUN returns a
#' list, or counts to be added). For a gen-backed object, a list with the
#' single result \code{FUN(as.matrix(x)[rows, ])}.
#' @param x A dartR object [required].
#' @param FUN Function of a genotype matrix (individuals x loci) [required].
#' @param rows Row indices, or NULL for all individuals [default NULL].
#' @return A list of per-block results, in locus order.
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.fbm.colblocks <- function(x, FUN, rows = NULL) {
  .gl_colblocks(x, FUN, rows = rows)
}

#' @name utils.fbm.rowwise
#' @title Adds a per-individual sum over loci in blocks of loci
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' Returns the same value as \code{FUN(as.matrix(x))}, where FUN is a sum
#' over loci per individual (e.g. \code{rowSums(is.na(m))}); FUN may return
#' a matrix with one row per individual (several sums). A file-backed object
#' is read in blocks of loci and the per-block results are added; sums of
#' counts are exact in double precision, so the values are identical. Only
#' for sums: for means of counts use \code{\link{utils.fbm.count.mean}} on
#' the summed counts.
#' @param x A dartR object [required].
#' @param FUN Function of a genotype matrix returning per-individual sums
#' [required].
#' @return A vector (or matrix) with one value (row) per individual.
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.fbm.rowwise <- function(x, FUN) {
  .gl_rowwise(x, FUN)
}

#' @name utils.fbm.count.mean
#' @title Means of counts, formed exactly as rowMeans() forms them
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' R's \code{mean()} of a logical vector, \code{rowMeans()} and
#' \code{colMeans()} divide a long-double sum by the count. For counts added
#' over blocks of loci, a plain \code{num / den} differs from that in the
#' last digit for about one count in 6,000 once \code{den} exceeds 2,048.
#' This returns, element by element, exactly the value \code{mean()} gives
#' for \code{num} TRUE among \code{den} values, so statistics counted in
#' blocks stay identical to those of the decoded matrix.
#' @param num Numerators: whole numbers (e.g. heterozygous loci per
#' individual) [required].
#' @param den Denominators: whole numbers, recycled to the length of num
#' (e.g. called loci per individual) [required].
#' @return A numeric vector with the names of num; NaN where den is 0.
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.fbm.count.mean <- function(num, den) {
  .gl_count_mean(num, den)
}

#' @name utils.fbm.reader
#' @title Reader of parts of the decoded genotype matrix
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' Returns a function \code{function(rows = NULL, cols = NULL)} giving
#' \code{as.matrix(x)[rows, cols, drop = FALSE]} (NULL: all rows or all
#' columns), with the dimnames and storage of \code{as.matrix()}. For a
#' file-backed object only the requested part is read from the backing
#' file; reading one individual touches the whole file (it is stored by
#' loci), so read few rows. For a gen-backed object the matrix is decoded
#' once and kept by the reader. Blocks of loci for a loop over the reader
#' come from \code{\link{utils.fbm.blocks}}.
#' @param x A dartR object [required].
#' @return A function of \code{rows} and \code{cols}.
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.fbm.reader <- function(x) {
  .gl_reader(x)
}

#' @name utils.fbm.blocks
#' @title Blocks of loci for reading a dartR object
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' A list of consecutive locus indices covering all loci, each block sized
#' by \code{options(bigstatsr.block.sizeGB)} for \code{nr} rows (the block
#' size dartR uses for file-backed objects).
#' @param x A dartR object [required].
#' @param nr Number of rows read per block [default nInd(x)].
#' @return A list of integer vectors.
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.fbm.blocks <- function(x, nr = nInd(x)) {
  .fbm_col_blocks(nLoc(x), nr)
}

#' @name utils.fbm.pair.counts
#' @title Pairwise genotype comparisons of all individuals of a file-backed
#' object
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' For each pair of individuals, the number of loci scored in both
#' (\code{L}) and, among those, the numbers that differ by one (\code{n1})
#' and by two allele copies (\code{n2}), counted over blocks of loci read
#' from the backing file (compiled with Rcpp at the first call). Genotypes
#' must be coded 0, 1, 2 or NA (diploid). The counts are whole numbers, so
#' they do not depend on the block size; distances built from them are
#' identical to those computed from the decoded matrix. Only for file-backed
#' objects; a gen-backed object gives an error.
#' @param x A file-backed dartR object [required].
#' @return A list of three integer vectors \code{L}, \code{n1}, \code{n2},
#' each in the order of a \code{dist} object (lower triangle by columns) of
#' length n(n - 1)/2 (12 bytes per pair).
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.fbm.pair.counts <- function(x) {
  .gl_pair_counts(x)
}
