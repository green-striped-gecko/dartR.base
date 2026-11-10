#' @name gl.print.history
#' @title Prints the history of a genlight object
#' @family environment
#'
#' @description
#' Prints the calls stored in the history of a genlight object
#' (\code{x@other$history}), one numbered entry per call, and returns them
#' as a table.
#'
#' @details
#' dartR functions that modify a genlight object add the call that produced
#' it to \code{x@other$history}. This function lists those calls in order.
#' Each entry keeps its position in the history, also when only some entries
#' are selected with \code{history}. Calls longer than 80 characters are
#' wrapped, with continuation lines indented by two spaces.
#'
#' To re-run the calls in a history, use \code{gl.play.history}.
#'
#' @param x A genlight object with a history. Not needed when
#' \code{history} is a history list [default NULL].
#' @param history Either a history list (such as \code{gl@other$history}),
#' or the numbers of the entries of \code{x@other$history} to print
#' (c(1, 3, 4) prints the first, third and fourth entries). If NULL, the
#' whole history of \code{x} is printed [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @return A data frame, returned invisibly, with one row per history entry:
#' \code{nr}, the position of the entry in the history, and \code{history},
#' the call as text. The table is printed at verbose 1 or higher.
#'
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' gl2 <- gl.filter.callrate(testset.gl, method = "loc", threshold = 0.9)
#' gl3 <- gl.filter.callrate(gl2, method = "ind", threshold = 0.95)
#' gl.print.history(gl3)
#' gl.print.history(gl3, history = c(1, 3))
#'
#' @export

gl.print.history <- function(x = NULL,
                             history = NULL,
                             verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname, verbose = verbose)

    # SCRIPT SPECIFIC CHECKS
    if (is.list(history)) {
        # A history list stands on its own; x is not needed
        hist2 <- history
        idx <- seq_along(hist2)
    } else {
        if (!is(x, "genlight")) {
            stop(error(
                "  Fatal Error: provide a genlight object as x, or a history",
                "list as history\n"
            ))
        }
        full <- x@other$history
        if (is.null(history)) {
            idx <- seq_along(full)
        } else {
            if (!is.numeric(history) || anyNA(history) ||
                any(history != round(history)) || any(history < 1) ||
                any(history > length(full))) {
                stop(error(
                    "  Fatal Error: history must be entry numbers from 1 to",
                    length(full), "\n"
                ))
            }
            idx <- history
        }
        hist2 <- full[idx]
    }

    # DO THE JOB

    # One line of text per entry
    calls <- vapply(hist2, function(h) {
        txt <- if (is.character(h)) h else deparse(h, width.cutoff = 500L)
        gsub("\\s+", " ", paste(txt, collapse = " "))
    }, character(1))
    dd <- data.frame(nr = idx, history = unname(calls),
                     stringsAsFactors = FALSE)

    if (nrow(dd) == 0) {
        if (verbose >= 2) {
            cat(warn("  Warning: no history entries found\n"))
        }
    } else if (verbose >= 1) {
        # Number, then the call; wrapped lines indented by two spaces
        w <- nchar(max(dd$nr))
        pad <- strrep(" ", w + 1)
        for (i in seq_len(nrow(dd))) {
            lines <- strwrap(dd$history[i], width = 80 - w - 1, exdent = 2)
            lead <- c(paste0(formatC(dd$nr[i], width = w), " "),
                      rep(pad, length(lines) - 1))
            cat(paste0(lead, lines, "\n"), sep = "")
        }
    }

    # FLAG SCRIPT END
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }

    invisible(dd)
}
