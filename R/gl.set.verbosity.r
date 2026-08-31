#' @name gl.set.verbosity
#' @title Sets the default verbosity level
#' @family environment
#' @description
#' dartR functions have a verbosity parameter that sets the level of reporting
#'  during the execution of the function. The verbosity level, set by parameter
#'  'verbose' can be one of verbose 0, silent or fatal errors; 1, begin and end;
#'  2, progress log; 3, progress and results summary; 5, full report. The
#'  default value for verbosity is stored in the r environment. This script sets
#'  the default value.
#' @param value Set the default verbosity to be this value: 0, silent only fatal
#' errors; 1, begin and end; 2, progress log; 3, progress and results summary;
#' 5, full report. An invalid value is coerced to 2 with a warning [default 2].
#' @return The verbosity value actually set, returned invisibly.
#' @export
#' @author Custodian: Arthur Georges (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl.set.verbosity(value = 2)

gl.set.verbosity <- function(value = 2) {
    # VALIDATE (before the banner, so an invalid value never reaches
    # the display machinery raw); invalid values coerce to the default
    # 2 with a warning, matching gl.check.verbosity
    if (is.null(value) || !is.numeric(value) || length(value) != 1 ||
        value < 0 || value > 5) {
        cat(
            warn(
                "Warning: Parameter value must be an integer in the range 0 to 5, set to 2\n"
            )
        )
        value <- 2
    }

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = value)

    # SET GLOBAL VERBOSITY
    options(dartR_verbose = value)
    if (value >= 2) {
        cat(report("  Global verbosity set to:", value, "\n"))
    }

    # FLAG SCRIPT END

    if (value >= 1) {
        cat(report("Completed:", funname, "\n"))
    }

    invisible(value)
}
