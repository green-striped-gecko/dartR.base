#' @name gl.propShared
#' @title Calculates a similarity matrix for individuals based on the
#' proportion of shared alleles
#' @family distance
#'
#' @description
#' This function calculates the proportion of shared alleles between each pair
#' of individuals and returns it as a similarity matrix.
#'
#' @details
#' At a locus, two diploid individuals with allele dosages g1 and g2 (0, 1 or
#' 2) share 2 - |g1 - g2| of their two alleles. The similarity is the mean
#' proportion of shared alleles over the loci called in both individuals,
#' 1 - mean(|g1 - g2|) / 2, so missing data are handled pair by pair.
#'
#' The similarity is the exact complement of
#' \code{\link{gl.dist.ind}} with \code{method = "manhattan"},
#' which this function uses for the calculation. Use
#' \code{as.dist(1 - gl.propShared(x))} to obtain a distance.
#'
#' Only SNP data are accepted: for SilicoDArT (presence/absence) data the
#' division by two does not apply.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @return A symmetric nInd(x) x nInd(x) similarity matrix with values
#' between 0 (no alleles shared at any locus) and 1 (identical at all loci
#' called in both individuals), and individual names as row and column names.
#' The diagonal is set to 1. Pairs of individuals with no locus called in
#' both get NA.
#'
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' res <- gl.propShared(bandicoot.gl)
#' res[1:5, 1:7] # show only a small part of the matrix
#'
#' @export

gl.propShared <- function(x,
                          verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  funname <- utils.flag.start(func = funname,
                   verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # DO THE JOB
  # 1 - Manhattan distance of gl.dist.ind is the proportion of shared alleles
  # in every dartR.base release; method "simple" matches only from PR #315
  res <- 1 - as.matrix(gl.dist.ind(x,
                                   method = "manhattan",
                                   type = "matrix",
                                   plot.display = FALSE,
                                   verbose = 0))
  diag(res) <- 1
  colnames(res) <- rownames(res) <- indNames(x)

  # FLAG SCRIPT END
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }

  return(res)
}
