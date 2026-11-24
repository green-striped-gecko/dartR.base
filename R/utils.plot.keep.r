#' @name utils.plot.keep
#' @title An internal function to keep the plots of the last dartR function
#' @family utilities
#'
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY
#' END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.
#'
#' Registers a plot (or another object that the function would save with
#' plot.file) so that \code{\link{gl.last.plot}} can return it. Call it in a
#' dartR function wherever a plot is made, whether or not it is displayed or
#' saved, with a name that is unique within the function.
#'
#' @param x The plot (usually a ggplot or patchwork object, or a base R
#' graphics plot recorded with utils.plot.record()) or other object. NULL
#' means that this call drew no plot: the plots of earlier functions are
#' then removed, so gl.last.plot() does not return an unrelated plot.
#' @param name Name of the object within the set kept for the calling
#' function, e.g. "plot", or "table" for a summary table [default "plot"].
#'
#' @details
#' All objects registered during one call of a function form one set, so a
#' function that makes several plots returns all of them. The set is
#' replaced when another function call registers a plot, so only the plots
#' of the most recent plotting function are kept. The calling function is
#' identified by a marker variable placed in its evaluation frame.
#'
#' A ggplot made inside a function keeps a reference to that function's
#' environment (its plot_env), which holds the function's arguments, so the
#' kept set holds them in memory until it is replaced or cleared with
#' gl.last.plot(clear = TRUE).
#'
#' @author Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#'
#' @keywords internal
#' @export
#' @return x, invisibly

utils.plot.keep <- function(x, name = "plot") {
  if (!(is.character(name) && length(name) == 1 && !is.na(name) &&
        nzchar(name))) {
    stop(error("Fatal Error: name must be a single non-empty character string.\n"))
  }
  caller <- parent.frame()
  # one marker per evaluation frame of the calling function: objects kept
  # during the same call share it and so form one set
  marker <- get0(".dartR_plot_call", envir = caller, inherits = FALSE)
  if (is.null(marker)) {
    .dartR_plots$counter <- .dartR_plots$counter + 1
    marker <- .dartR_plots$counter
    assign(".dartR_plot_call", marker, envir = caller)
  }
  if (!identical(.dartR_plots$marker, marker)) {
    call <- sys.call(-1)
    fun <- if (is.null(call)) "" else paste(deparse(call[[1]]), collapse = "")
    .dartR_plots$marker <- marker
    .dartR_plots$fun <- fun
    .dartR_plots$objects <- list()
  }
  # NULL (no plot drawn) leaves the new set empty
  if (!is.null(x)) {
    .dartR_plots$objects[[name]] <- x
  }
  invisible(x)
}

#' @name utils.plot.record
#' @title An internal function to record a plot drawn with base R graphics
#' @family utilities
#'
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY
#' END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.
#'
#' Records the plot on the current graphics device with
#' grDevices::recordPlot(), so that a function that draws with base R
#' graphics (e.g. a heatmap or a tree) can keep it with
#' \code{\link{utils.plot.keep}}. Call it right after the plot is drawn
#' completely.
#'
#' @details
#' A recorded plot is redrawn with print() or grDevices::replayPlot(); it
#' cannot be changed like a ggplot. Devices that write to a file (pdf, png,
#' as in Rscript) keep no display list, so nothing can be recorded there and
#' NULL is returned; the screen devices of R and RStudio can be recorded.
#'
#' @author Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#'
#' @keywords internal
#' @export
#' @return A recordedplot object, or NULL if no device is open or the
#' device keeps no display list.

utils.plot.record <- function() {
  if (grDevices::dev.cur() == 1) {
    return(NULL)
  }
  rp <- tryCatch(suppressWarnings(grDevices::recordPlot()),
                 error = function(e) NULL)
  # an empty display list: the device does not record (file devices)
  if (is.null(rp) || length(rp[[1]]) == 0) {
    return(NULL)
  }
  rp
}

# The plots kept by utils.plot.keep(); read by gl.last.plot()
.dartR_plots <- new.env(parent = emptyenv())
.dartR_plots$counter <- 0
.dartR_plots$marker <- NULL
.dartR_plots$fun <- ""
.dartR_plots$objects <- list()
