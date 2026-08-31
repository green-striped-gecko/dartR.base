#' @name gl.He
#' @title Estimates expected heterozygosity for each locus
#' @family utilities

#' @description
#' Calculates the expected heterozygosity at each locus under Hardy-Weinberg
#' proportions, He = 2p(1-p), with the allele frequency p estimated from the
#' mean allele dosage pooled across all individuals in the genlight object.

#' @details
#' This is a silent pure accessor: it returns the per-locus vector directly
#' (visibly), produces no console output and does not modify the input. The
#' estimate carries no sample-size correction — for unbiased expected
#' heterozygosity (uHe) averaged within populations, use
#' \code{\link{gl.report.heterozygosity}}.
#'
#' Loci with all genotypes missing return NaN. The function applies to SNP
#' genotype data only.

#' @param gl A genlight object containing SNP genotypes [required].

#' @return A named vector with the expected heterozygosity for each locus.

#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})

#' @examples
#' pp <- possums.gl[1:30,1:20]
#' if (isTRUE(getOption("dartR_fbm"))) pp <- gl.gen2fbm(pp)
#' gl.He(pp)

#' @seealso \code{\link{gl.Ho}}, \code{\link{gl.alf}},
#' \code{\link{gl.report.heterozygosity}}
#' @export

gl.He <- function(gl) {
  utils.check.datatype(gl, accept = "SNP", verbose = 0)
  alf <- colMeans(as.matrix(gl), na.rm = T) / 2
  out <- alf * (1 - alf) * 2
  return(out)
}
