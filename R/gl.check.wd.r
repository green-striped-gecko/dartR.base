#' @name gl.check.wd
#' @title Checks the global working directory
#' @family environment
#'
#' @description
#' The working directory can be set in one of two ways -- (a) explicitly by the
#' user by passing a value using the parameter plot.dir in a function, or (b) by
#' setting the working directory globally as part of the r environment
#' (gl.set.wd). In accordance with CRAN policy, the default is tempdir().

#' @param wd path to the working directory [default NULL; resolves to the
#' dartR_wd option if set with gl.set.wd(), otherwise to tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @examples
#' gl.check.wd()
#'
#' @author Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#'
#' @export
#' @return the working directory

gl.check.wd <- function(
    wd = NULL,
    verbose=NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)
  # DO THE JOB
  # SET wd or GET it from global
  if (is.null(wd)) {
    # If wd is not provided, check if it's set in options()
    if (is.null(options()$dartR_wd)) {
      # If not set in options(), set wd to tempdir()
      wd <- tempdir()
    } else {
      # If set in options(), use that value for wd
      wd <- options()$dartR_wd
    }
  } else {
    # If wd is provided: it must be a single existing directory path. Use &&
    # so dir.exists() is not evaluated on a non-character wd, and require
    # length 1 so the if() condition is scalar (a non-character, NA, empty or
    # multi-element wd falls through to the tempdir fallback rather than
    # raising an opaque error).
    if (is.character(wd) && length(wd) == 1 && !is.na(wd) && dir.exists(wd)) {
      # a valid directory path: keep wd as supplied
    } else {
      # not a valid directory path: warn (gated) and fall back to tempdir
      if (verbose >= 1) {
        cat(
          warn(
            "Warning: The path to the working directory does not exist! Set to tempdir().\n"
          )
        )
      }
      wd <- tempdir()
    }
  }
  if(verbose >= 2){cat(report("  Working directory:",wd,"\n"))}
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

    return(wd)

}
