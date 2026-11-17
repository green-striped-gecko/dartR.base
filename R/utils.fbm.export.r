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

#' @name utils.fbm.grm
#' @title Additive relationship matrix of a file-backed object, as
#' rrBLUP::A.mat with its default settings, in blocks of loci
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' Computes the realized additive relationship matrix of
#' \code{rrBLUP::A.mat(as.matrix(x) - 1, min.MAF, max.missing)} with A.mat's
#' default mean imputation and no shrinkage (Endelman & Jannink 2012;
#' VanRaden 2008), reading the genotypes in two passes over blocks of loci,
#' so the genotype matrix is never decoded as a whole.
#'
#' As in A.mat, with genotypes g (0, 1, 2) and \code{X = g - 1}: per locus
#' the frequency p is \code{mean(g, na.rm = TRUE) / 2} (computed per column
#' exactly as A.mat computes it) and the missing fraction is the share of
#' NA. Loci with \code{min(p, 1 - p) >= min.MAF} and missing fraction
#' \code{<= max.missing} are used (m loci). \code{W = g - 2p}, with missing
#' genotypes set to 0 (mean imputation) and loci fixed in the sample set to
#' 0, and \code{A = tcrossprod(W) / (2 * mean(p * (1 - p))) / m}, where the
#' mean is over the m loci.
#'
#' The frequencies, filter and scaling are identical to A.mat; the
#' cross-product is added block by block, so A agrees with A.mat to about
#' 1e-13 (relative) rather than bit for bit. Memory: the n x n result and
#' one block of loci.
#' @param x A file-backed dartR object (diploid SNP) [required].
#' @param min.MAF Minimum minor allele frequency [default 1/(2n), as A.mat].
#' @param max.missing Maximum missing fraction of a locus [default
#' 1 - 1/(2n), as A.mat].
#' @return An n x n matrix with indNames(x) as row and column names.
#' @references Endelman, J.B. and Jannink, J.-L. (2012) Shrinkage estimation
#' of the realized relationship matrix. G3 2:1405-1413.
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.fbm.grm <- function(x, min.MAF = NULL, max.missing = NULL) {
  if (!.has_fbm(x)) {
    stop(error("  Fatal Error: utils.fbm.grm() needs a file-backed object;",
               "use rrBLUP::A.mat for a gen-backed one\n"))
  }
  n <- nInd(x)
  blocks <- .fbm_col_blocks(nLoc(x), n)
  rd <- .gl_reader(x)
  # Pass 1: per locus the number of NA and the mean genotype, as A.mat:
  # apply(X + 1, 2, mean) on doubles (two-pass mean() of a double vector)
  s1 <- do.call(cbind, lapply(blocks, function(cols) {
    g <- rd(cols = cols)
    storage.mode(g) <- "double"
    rbind(colSums(is.na(g)),
          apply(g, 2, function(v) mean(v, na.rm = TRUE)))
  }))
  frac.missing <- s1[1, ] / n
  missing <- max(frac.missing) > 0
  freq <- s1[2, ] / 2
  MAF <- pmin(freq, 1 - freq)
  if (is.null(min.MAF)) min.MAF <- 1 / (2 * n)
  if (is.null(max.missing)) max.missing <- 1 - 1 / (2 * n)
  use <- (MAF >= min.MAF) & (frac.missing <= max.missing)
  use[is.na(use)] <- FALSE
  markers <- which(use)
  m <- length(markers)
  var.A <- 2 * mean(freq[markers] * (1 - freq[markers]))
  # loci fixed in the sample: A.mat replaces their genotypes (also NA) by
  # 2p - 1, so their column of W is 0
  mono <- freq * (1 - freq) == 0
  # Pass 2: W W' added over blocks of the loci used
  A <- matrix(0, n, n)
  for (cols in blocks) {
    k <- cols[use[cols]]
    if (length(k) == 0) next
    g <- rd(cols = k)
    storage.mode(g) <- "double"
    W <- sweep((g - 1) + 1, 2, 2 * freq[k])
    if (missing) W[is.na(W)] <- 0
    W[, mono[k]] <- 0
    A <- A + tcrossprod(W)
  }
  A <- A / var.A / m
  dimnames(A) <- list(indNames(x), indNames(x))
  A
}

#' @name utils.fbm.amat
#' @title Additive relationship matrix by rrBLUP::A.mat for gen- and
#' file-backed objects
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' Computes \code{rrBLUP::A.mat(as.matrix(x) - 1, ...)}, the matrix used by
#' gl.grm (dartR.captive) and gl.grm2 (dartR.spatial):
#' \itemize{
#' \item gen-backed object: rrBLUP::A.mat, unchanged;
#' \item file-backed object with A.mat's default method (mean imputation,
#' no shrinkage, no imputed matrix returned; \code{min.MAF},
#' \code{max.missing}, \code{n.core} and \code{tol} may be given):
#' \code{\link{utils.fbm.grm}}, in blocks of loci without decoding the
#' genotype matrix (agrees with A.mat to about 1e-13);
#' \item file-backed object with any other A.mat option (\code{shrink},
#' \code{impute.method = "EM"}, \code{return.imputed}): these need the whole
#' matrix, so the object is converted to a genlight (\code{gl.fbm2gen}) with
#' a warning that states the memory A.mat will use, and A.mat is run (it
#' stops instead if that exceeds the limit of
#' \code{\link{utils.fbm.decode}}).
#' }
#' @param x A dartR object (diploid SNP) [required].
#' @param dots A named list of arguments for rrBLUP::A.mat [default list()].
#' @param verbose Verbosity of the warning about converting [default 2].
#' @return The matrix returned by rrBLUP::A.mat (or a list, with
#' \code{return.imputed = TRUE}).
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.fbm.amat <- function(x, dots = list(), verbose = 2) {
  if (!requireNamespace("rrBLUP", quietly = TRUE)) {
    stop(error("  Package rrBLUP needed for this function to work. Please",
               "install it.\n"))
  }
  if (!.has_fbm(x)) {
    return(do.call(rrBLUP::A.mat, c(list(as.matrix(x) - 1), dots)))
  }
  nm <- names(dots)
  if (is.null(nm)) nm <- rep("", length(dots))
  default <- all(nm %in% c("min.MAF", "max.missing", "n.core", "tol",
                           "impute.method", "shrink", "return.imputed")) &&
    (is.null(dots$impute.method) ||
       identical(toupper(dots$impute.method), "MEAN")) &&
    (is.null(dots$shrink) || isFALSE(dots$shrink)) &&
    (is.null(dots$return.imputed) || isFALSE(dots$return.imputed))
  if (default) {
    return(utils.fbm.grm(x, min.MAF = dots$min.MAF,
                         max.missing = dots$max.missing))
  }
  # A.mat holds about four n x p matrices of doubles (input, filtered and
  # centred copies, frequency matrix) besides the n x n result
  gb <- 4 * 8 * as.numeric(nInd(x)) * nLoc(x) / 1024^3
  # stop before exhausting memory (limit: option dartR_fbm_max_decode_gb)
  .fbm_decode_guard(x, what = "rrBLUP::A.mat with these options",
                    bytes = 8, copies = 4, verbose = 0)
  if (verbose >= 1) {
    cat(warn("  Warning: the rrBLUP::A.mat option(s)",
             paste(setdiff(nm, c("min.MAF", "max.missing", "n.core", "tol")),
                   collapse = ", "),
             "need the whole genotype matrix; converting the file-backed",
             "object to a genlight. A.mat will use about",
             signif(gb, 2), "Gb of memory\n"))
  }
  x <- gl.fbm2gen(x, verbose = 0)
  do.call(rrBLUP::A.mat, c(list(as.matrix(x) - 1), dots))
}

#' @name utils.fbm.decode
#' @title Decodes the whole genotype matrix, with a memory guard for
#' file-backed objects
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' For functions whose algorithm needs the whole genotype matrix in memory.
#' Returns \code{as.matrix(x)}. For a gen-backed object nothing else
#' happens. For a file-backed object the memory the matrix needs is
#' estimated first (\code{nInd(x) * nLoc(x) * bytes}): it is reported at
#' verbose >= 2, and if it exceeds
#' \code{getOption("dartR_fbm_max_decode_gb", 4)} Gb the function stops
#' with an informative error instead of exhausting memory. The limit can
#' be raised, e.g. \code{options(dartR_fbm_max_decode_gb = 16)}.
#' @param x A dartR object [required].
#' @param what Name of the calling function, for the messages
#' [default "this function"].
#' @param bytes Bytes per genotype of the matrix the caller works with
#' (4 for integers, 8 for doubles) [default 8].
#' @param copies Number of such matrices the caller holds at once
#' [default 1].
#' @param verbose Verbosity [default 2].
#' @return The decoded genotype matrix, as \code{as.matrix(x)}.
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.fbm.decode <- function(x, what = "this function", bytes = 8,
                             copies = 1, verbose = 2) {
  if (!.has_fbm(x)) return(as.matrix(x))
  .fbm_decode_guard(x, what = what, bytes = bytes, copies = copies,
                    verbose = verbose)
  as.matrix(x)
}

## Estimates the memory of decoding a file-backed x (copies matrices of
## 'bytes' per genotype), reports it at verbose >= 2 and stops above
## getOption("dartR_fbm_max_decode_gb", 4) Gb
.fbm_decode_guard <- function(x, what, bytes = 8, copies = 1, verbose = 2) {
  gb <- copies * bytes * as.numeric(nInd(x)) * nLoc(x) / 1024^3
  limit <- getOption("dartR_fbm_max_decode_gb", 4)
  if (!is.numeric(limit) || length(limit) != 1 || is.na(limit)) {
    stop(error("  Fatal Error: option dartR_fbm_max_decode_gb must be a",
               "single number (Gb)\n"))
  }
  if (gb > limit) {
    stop(error("  Fatal Error:", what, "needs the whole genotype matrix in",
               "memory: about", signif(gb, 2), "Gb for", nInd(x),
               "individuals and", nLoc(x), "loci, more than the limit of",
               limit, "Gb (option dartR_fbm_max_decode_gb). Raise the",
               "limit if enough memory is available, e.g.",
               paste0("options(dartR_fbm_max_decode_gb = ",
                      ceiling(gb * 1.5), "),"),
               "or reduce the data (filter loci or individuals)\n"))
  }
  if (verbose >= 2) {
    cat(report("  ", what, "needs the whole genotype matrix in memory:",
               "decoding the file-backed object (about", signif(gb, 2),
               "Gb)\n"))
  }
  invisible(gb)
}
