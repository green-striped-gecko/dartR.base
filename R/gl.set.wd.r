#' @name gl.set.wd
#' @title Sets the default working directory
#' @family environment

#' @description
#' Many dartR functions have a plot.dir parameter which is used to save output to (e.g. ggplots as rds files)
#'  With this functions users can set the working directory globally so it is used in all functions, without setting is explicitely. 
#'  The value for wd is stored in the r environment and if not set defaults to tempdir(). This script sets the default value.
#' @param wd Set the path to the wd directory globally to be used by all functions if not set explicitely in the function.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].

#' @author Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#' #set to current working directory
#' wd <- gl.set.wd(wd=getwd())
#' 
#' @export
#' @return path the the working directory [set for all functions]

gl.set.wd <- function(
    wd = tempdir(),
    verbose=NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbose = verbose)
  # DO THE JOB
  if (!is.null(wd) &    # Check if 'wd' argument is not NULL
      is.character(wd) &    # Check if 'wd' is a character string
      dir.exists(wd)){    # Check if the directory specified by 'wd' exists
    
    options(dartR_wd = wd)    # Set the 'dartR_wd' option to the specified directory
  }
  if(verbose >= 2){cat(report("  Global working directory set to",wd,"\n"))}
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
    return(wd)
}
