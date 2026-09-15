#' @name gl.set.wd
#' @title Sets the default working directory
#' @family environment

#' @description
#' Many dartR functions have a plot.dir parameter which is used to save output
#' to (e.g. ggplots as rds files). With this function users can set the working
#' directory globally so it is used in all functions, without setting it
#' explicitly. The value for wd is stored in the r environment and if not set
#' defaults to tempdir(). This function sets the default value.

#' @details
#' The nominated directory must already exist; gl.set.wd() does not create it.
#' If the path does not exist (or is not a single directory path) the function
#' stops with an error and the global working directory is left unchanged.

#' @param wd Path to the directory to set globally as the working directory,
#' to be used by all functions when not set explicitly in the function
#' [default tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @author Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' #set to current working directory
#' wd <- gl.set.wd(wd=getwd())
#'
#' @export
#' @return path to the working directory [set for all functions]

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

  # FUNCTION SPECIFIC ERROR CHECKING
  # wd must be a single existing directory path. Use && so dir.exists() is
  # not evaluated on a non-character/NULL wd, and require length 1 so the
  # condition is scalar. Fail loudly rather than silently reporting success
  # when nothing would be set.
  if (!(is.character(wd) && length(wd) == 1 && !is.na(wd) && dir.exists(wd))) {
    stop(error(
      "Fatal Error: the nominated working directory does not exist (or is not",
      " a single directory path); the global working directory was not",
      " changed. Create the directory first, or pass an existing path.\n"
    ))
  }

  # DO THE JOB
  options(dartR_wd = wd)    # Set the 'dartR_wd' option to the specified directory

  if(verbose >= 2){cat(report("  Global working directory set to",wd,"\n"))}
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

    return(wd)
}
