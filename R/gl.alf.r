#' @name gl.alf
#' @title Calculates the frequency of the reference and alternate allele for
#' each locus
#' @family utilities

#' @description
#' Calculates the frequency of the reference allele (alf1) and the alternate
#' allele (alf2) at each locus, pooled across all individuals in the genlight
#' object.

#' @details
#' This is a lightweight accessor: it returns the per-locus frequencies
#' directly (visibly), produces no console output and does not modify the
#' input. The alternate allele frequency is the mean genotype score at the
#' locus divided by 2, computed with na.rm = TRUE, and alf1 is its
#' complement, so the two columns sum to 1 at every locus with at least one
#' scored genotype. A locus with no scored genotypes returns NaN in both
#' columns.
#'
#' The frequencies are unrounded, and they are computed in a single pass over
#' the genotype matrix without splitting the object by population, which
#' makes this the fast path used inside loops and by \code{utils.recalc.maf}.
#' \code{\link{gl.allele.freq}} with \code{simple = TRUE} returns the same
#' quantity rounded to 4 decimal places, and additionally offers the
#' per-population and population x locus aggregations that gl.alf does not.
#'
#' Row names are the locus names, in locus order. Duplicate locus names are
#' disambiguated with \code{make.unique()} so that the row keys always
#' correspond positionally to the loci of the input.
#'
#' The function applies to SNP genotype data only. Tag presence/absence
#' (SilicoDArT) data has ploidy 1 and no alternate allele, so it is rejected
#' rather than divided by the SNP ploidy.

#' @param x Name of the genlight object [required].

#' @return A data.frame with one row per locus, in locus order, and two
#' columns: alf1, the frequency of the reference allele, and alf2, the
#' frequency of the alternate allele. Loci with all genotypes missing return
#' NaN in both columns.

#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#'
#' #test fbm
#' if (isTRUE(getOption("dartR_fbm"))) possums.gl <- gl.gen2fbm(possums.gl)
#' #for the first 10 loci only
#' gl.alf(possums.gl[,1:10])
#' barplot(t(as.matrix(gl.alf(possums.gl[,1:10]))))
#' gl.allele.freq(possums.gl[,1:10],simple=TRUE)
#' barplot(t(as.matrix(gl.allele.freq(possums.gl[,1:10],simple=TRUE))))

#' @seealso \code{\link{gl.allele.freq}}, \code{\link{gl.Ho}},
#' \code{\link{gl.He}}
#' @export

gl.alf <- function(x) {
  utils.check.datatype(x, accept = "SNP", verbose = 0)
  alf <- colMeans(as.matrix(x), na.rm = TRUE) / 2
  # data.frame() accepts the colMeans names as row names only when they are
  # unique, and otherwise substitutes 1:n without warning. Callers read
  # rownames() as locus names, so set them explicitly from locNames() and
  # keep positional correspondence with the loci of the input.
  out <- data.frame(alf1 = 1 - alf, alf2 = alf, row.names = NULL)
  loc.names <- locNames(x)
  if (!is.null(loc.names)) {
    rownames(out) <- make.unique(loc.names)
  }
  return(out)
}
