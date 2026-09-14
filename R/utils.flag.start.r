#' @name utils.flag.start
#' @title A utility script to flag the start of a script
#' @family utilities

#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.

#' @param func Name of the function that is starting [required].
#' @param build Name of the build [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity]
#'  

#' @author Custodian: Arthur Georges -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#'  
#' @keywords internal
#' @export
#' @return The calling function name, invisibly
 
# Build Version v.2023.2

utils.flag.start <- function(func = NULL,
                             build = NULL,
                             verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    if (is.null(func)) {
        stop(error("Fatal Error: The calling function must be specified.\n"))
    }
    if (verbose >= 1) {
        if (verbose == 5 && !is.null(build)) {
            cat(
                report(
                    "Starting",
                    func,
                    "\n[dartR.base vers.",
                    packageVersion("dartR.base"),
                    "Build =",
                    build,
                    "]\n"
                )
            )
        } else {
            cat(report("Starting", func, "\n"))
        }
    }
    invisible(func)
}
