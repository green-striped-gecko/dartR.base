#' Estimates observed Heterozygosity

#' @param gl A genlight object [required]
#' @return A simple vector whit Ho for each loci
#' @export
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples 
#' pp <- possums.gl[1:30,1:20]
#' if (isTRUE(getOption("dartR_fbm"))) pp <- gl.gen2fbm(pp)
#' gl.Ho(pp)

gl.Ho <- function(gl) {
  out <- colMeans(as.matrix(gl) == 1, na.rm = T)
  return(out)
}
