#' Estimates expected Heterozygosity

#' @param gl A genlight object [required]
#' @return A simple vector whit Ho for each loci
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples 
#' pp <- possums.gl[1:30,1:20]
#' if (isTRUE(getOption("dartR_fbm"))) pp <- gl.gen2fbm(pp)
#' gl.He(pp)

gl.He <- function(gl) {
  alf <- colMeans(as.matrix(gl), na.rm = T) / 2
  out <- alf * (1 - alf) * 2
  return(out)
}
