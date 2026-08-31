#' @name gl.Ho
#' @title Estimates observed heterozygosity for each locus
#' @family utilities

#' @description
#' Calculates the observed heterozygosity at each locus, pooled across all
#' individuals in the genlight object, as the proportion of individuals
#' scored as heterozygous at that locus.

#' @details
#' This is a silent pure accessor: it returns the per-locus vector directly
#' (visibly), produces no console output and does not modify the input. For
#' heterozygosities averaged within populations or by individual, with
#' sample-size corrections and plots, use
#' \code{\link{gl.report.heterozygosity}}.
#'
#' Loci with all genotypes missing return NaN. The function applies to SNP
#' genotype data only.

#' @param gl A genlight object containing SNP genotypes [required].

#' @return A named vector with the observed heterozygosity for each locus.

#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber (bugs? Post to
#' \url{https://groups.google.com/d/forum/dartr})

#' @examples
#' pp <- possums.gl[1:30,1:20]
#' if (isTRUE(getOption("dartR_fbm"))) pp <- gl.gen2fbm(pp)
#' gl.Ho(pp)

#' @seealso \code{\link{gl.He}}, \code{\link{gl.alf}},
#' \code{\link{gl.report.heterozygosity}}
#' @export

gl.Ho <- function(gl) {
  utils.check.datatype(gl, accept = "SNP", verbose = 0)
  out <- colMeans(as.matrix(gl) == 1, na.rm = T)
  return(out)
}
