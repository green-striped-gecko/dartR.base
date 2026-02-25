#' @name gl.alf
#' @title
#' Calculates allele frequency of the first and second allele for each locus
#' A very simple function to report allele frequencies
#' @family utilities

#' @param x Name of the genlight object [required].
#' @author Bernd Gruber (bugs? Post to
#'  \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' 
#' #test fbm
#' if (isTRUE(getOption("dartR_fbm"))) possums.gl <- gl.gen2fbm(possums.gl)
#' #for the first 10 loci only
#' gl.alf(possums.gl[,1:10])
#' barplot(t(as.matrix(gl.alf(possums.gl[,1:10]))))
#' gl.allele.freq(possums.gl[,1:10],simple=TRUE)
#' barplot(t(as.matrix(gl.allele.freq(possums.gl[,1:10],simple=TRUE))))
#' @export
#' @return A simple data.frame with ref (reference allele), alt (alternate allele).

gl.alf <- function(x) {
  #cat(warn("Deprecated: Please use gl.allele.freq(x,simple=TRUE)\n"))
  alf <- colMeans(as.matrix(x), na.rm = TRUE) / 2
  out <- data.frame(alf1 = 1 - alf, alf2 = alf)
  return(out)
}
