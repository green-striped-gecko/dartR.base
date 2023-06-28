#' @name gl.setwd
#' @title Sets the default working directory
#' @description
#' Many dartR functions have a plot.dir parameter which is used to save output to (e.g. ggplots as rds files)
#'  With this functions users can set the working directory globally so it is used in all functions, without setting is explicitely. 
#'  The value for wd is stored in the r environment and if not set defaults to tempdir(). This script sets the default value.
#' @param wd Set the path to the wd directory globally to be used by all functions if not set explicitely in the function.
#' @return path the the working directory [set for all functions]
#' @export
#' @author Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' #set to current working directory
#' gl <- gl.setwd(wd=getwd())

gl.setwd <- function(wd = tempdir()) {
    # SET GLOBAL VERBOSITY
    if (!is.null(wd) &
        is.character(wd) & dir.exists(wd)){
        options(dartR_wd = wd)
    }
    return(wd)
}
