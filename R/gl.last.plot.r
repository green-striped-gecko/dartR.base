#' @name gl.last.plot
#' @title Returns the plot(s) of the last dartR function that made a plot
#' @family graphics
#'
#' @description
#' Many dartR functions draw a plot but return something else, for example
#' the filtered genlight object. gl.last.plot() returns the plot that the
#' most recent such function made, so it can be changed or saved, without
#' saving it to a file first (see the plot.file and plot.dir arguments of
#' those functions).
#'
#' @param which Name or number of one object to return when the function
#' made several (see Details). If NULL, all are returned [default NULL].
#' @param clear If TRUE, the kept plots are removed from memory and NULL is
#' returned [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @details
#' The plot is kept whether or not it was displayed (plot.display) or saved
#' (plot.file). Only the plots of the most recent plotting function are
#' kept; the next function that makes a plot replaces them.
#'
#' \strong{One plot.} The plot itself is returned, usually a ggplot. Several
#' panels drawn as one figure (e.g. a histogram above a boxplot) form one
#' patchwork object: change one panel with p[[1]], all panels with the &
#' operator (see the examples and \code{\link{gl.plot.help}}).
#'
#' \strong{Several objects.} Some functions make several plots, or a plot
#' and a summary table (the objects that plot.file would save in separate
#' files). They are returned as a named list, invisibly; the names are
#' reported at verbose >= 2. Use which to get one of them.
#'
#' A kept plot holds the data it was drawn from in memory until it is
#' replaced; use gl.last.plot(clear = TRUE) to free it.
#'
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' t1 <- gl.filter.callrate(testset.gl, threshold = 0.95, verbose = 0)
#' p <- gl.last.plot()
#' # the call rate plots are two panels: change the first one only
#' p[[1]] <- p[[1]] + ggplot2::labs(title = "Call rate before filtering")
#' # or the text size of both panels
#' p & ggplot2::theme(text = ggplot2::element_text(size = 14))
#' # code to change axes, colours, symbols or single points
#' gl.plot.help(p)
#'
#' @seealso \code{\link{gl.plot.help}}, \code{\link{theme_dartR}}
#' @export
#' @return The plot, or a named list of the objects when the function made
#' several (invisibly), or NULL when clear = TRUE.

gl.last.plot <- function(which = NULL,
                         clear = FALSE,
                         verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  funname <- utils.flag.start(func = funname,
                              build = "v.2026.1",
                              verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (!(is.logical(clear) && length(clear) == 1 && !is.na(clear))) {
    stop(error("Fatal Error: clear must be TRUE or FALSE.\n"))
  }

  # DO THE JOB
  if (clear) {
    .dartR_plots$marker <- NULL
    .dartR_plots$fun <- ""
    .dartR_plots$objects <- list()
    if (verbose >= 2) {
      cat(report("  Kept plots removed from memory\n"))
    }
    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }
    return(invisible(NULL))
  }

  obj <- .dartR_plots$objects
  if (length(obj) == 0) {
    stop(error("Fatal Error: no dartR function has made a plot yet in this",
               "session (or the plots were cleared).\n"))
  }
  fun <- .dartR_plots$fun
  nms <- names(obj)

  if (!is.null(which)) {
    ok <- (is.character(which) && length(which) == 1 && which %in% nms) ||
      (is.numeric(which) && length(which) == 1 && !is.na(which) &&
         which == round(which) && which >= 1 && which <= length(obj))
    if (!ok) {
      stop(error("Fatal Error: which must be one of",
                 paste0("'", nms, "'", collapse = ", "),
                 "or a number from 1 to", length(obj), "\n"))
    }
    out <- obj[[which]]
    if (verbose >= 2) {
      cat(report("  Returning", paste0("'", if (is.numeric(which)) nms[which] else which, "'"),
                 "made by", paste0(fun, "()"), "\n"))
    }
    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }
    return(out)
  }

  if (length(obj) == 1) {
    if (verbose >= 2) {
      cat(report("  Returning the plot made by", paste0(fun, "()"), "\n"))
    }
    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }
    return(obj[[1]])
  }

  if (verbose >= 2) {
    cat(report(paste0("  ", fun, "()"), "made", length(obj), "objects:",
               paste0("$", nms, collapse = ", "),
               "\n  Returning them as a named list; use which to get one of them\n"))
  }
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  invisible(obj)
}
