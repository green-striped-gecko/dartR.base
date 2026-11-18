#' @name gl.check.verbosity
#' @title Returns the verbosity level a function should use
#' @family environment
#'
#' @description
#' Resolves the verbosity level for a dartR function. Most dartR functions
#' call it first, as \code{verbose <- gl.check.verbosity(verbose)}.
#'
#' @details
#' The verbosity comes from one of three sources, in this order of
#' precedence:
#' \enumerate{
#'   \item the value passed by the user in the function's \code{verbose}
#'   argument;
#'   \item the global default set with \code{\link{gl.set.verbosity}},
#'   stored as the R option \code{dartR_verbose};
#'   \item the value 2, when neither of the above is set.
#' }
#' A valid value is a single number from 0 to 5. Any other value (out of
#' range, non-numeric, \code{NA}, or a vector of length other than one),
#' whether passed as \code{x} or found in the global option, prints a
#' warning and is replaced by 2.
#'
#' @param x User requested level of verbosity, a single number from 0 to 5,
#' or NULL to use the global default [default NULL].
#'
#' @return The verbosity level to use, a single number from 0 to 5.
#'
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' gl.check.verbosity()
#' gl.check.verbosity(3)
#'
#' @export

# Version v.2023.3

gl.check.verbosity <- function(x = NULL) {
    # SET VERBOSITY or GET it from global; an explicit value wins
    source <- "verbose"
    if (is.null(x)) {
        x <- getOption("dartR_verbose")
        source <- "option dartR_verbose"
        if (is.null(x)) {
            return(2)
        }
    }

    # A computed value can be NA or a vector, so check length and NA
    # before the range test; otherwise `if` stops with an R-internal error
    if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 0 || x > 5) {
        cat(warn(paste0(
            "Warning: ", source, " must be a single number from 0 to 5",
            " (received ", paste(deparse(x), collapse = ""), "); set to 2\n"
        )))
        return(2)
    }

    return(x)
}
