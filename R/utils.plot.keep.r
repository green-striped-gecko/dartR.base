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
#' A ggplot made inside a function refers to that function's environment
#' (its plot_env and the environments of its aes() mappings), which holds
#' all objects of the function, e.g. the genlight, a decoded genotype matrix
#' or a filtered copy. To keep these from staying in memory, the kept copy
#' of a ggplot (also each panel of a patchwork) is "slimmed": every
#' reference to the function's environment is replaced by a small
#' environment that holds only the objects the plot refers to by name. The
#' plot is unchanged; the object returned by the function is not modified.
#' If slimming fails, the plot is kept as it is.
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
    # keep a slimmed copy: it does not hold the function's objects in memory
    .dartR_plots$objects[[name]] <- utils.plot.slim(x, caller)
  }
  invisible(x)
}

# ---- slimming of kept plots (internal, not exported) ------------------------

# A copy of a ggplot or patchwork object that no longer refers to the
# evaluation frame of the function that made it (frame). Every reference
# to frame (or an environment inside it) is replaced by one small
# environment holding only the objects the plot refers to by name, with
# the parent of frame (usually the package namespace) as its parent, so
# functions still resolve. Layers, facet and scales are copied before they
# are changed (they are environments shared with the original plot), so the
# plot returned by the function is not modified. Any other object, e.g. a
# table or a recorded base R plot, is returned as it is. If anything fails,
# x is returned unchanged.
utils.plot.slim <- function(x, frame) {
  if (!(inherits(x, "ggplot") || inherits(x, "patchwork"))) {
    return(x)
  }
  # if slimming fails, the plot is kept unslimmed; the reason (message and
  # the calls that led to it) is stored in .dartR_plots$slim_error, so a
  # failure can be found (e.g. by a test) instead of passing unnoticed
  where <- NULL
  out <- tryCatch(
    withCallingHandlers(utils.plot.slim.do(x, frame), error = function(e) {
      # the calls up to the error, without the two of this handler
      calls <- utils::head(sys.calls(), -2)
      where <<- vapply(utils::tail(calls, 8), function(cl) {
        substr(paste(deparse(cl), collapse = " "), 1, 120)
      }, character(1))
    }),
    error = function(e) {
      .dartR_plots$slim_error <- c(paste("Error:", conditionMessage(e)), where)
      x
    })
  out
}

utils.plot.slim.do <- function(x, frame) {
  # The slimmed plot hangs from the first named environment above frame (for
  # a dartR function its namespace, otherwise e.g. the global environment).
  # Unnamed environments in between (e.g. that of a function defined inside
  # another function, or in a test) would be saved with the plot together
  # with everything in them; the objects the plot needs from them are copied
  # into slim below.
  anchor <- parent.env(frame)
  while (!nzchar(environmentName(anchor)) && !identical(anchor, emptyenv())) {
    anchor <- parent.env(anchor)
  }
  slim <- new.env(parent = anchor)
  # the namespace (or global environment) of the function that keeps the plot
  home <- environmentName(anchor)
  # Is e the evaluation frame of a function (or an environment inside one)
  # whose code lives in a dartR package or in the keeping function's home?
  # This covers frame itself and the frames of dartR functions called by it
  # that made parts of the plot (e.g. gl.smearplot() panels). Frames of other
  # packages' functions (e.g. label closures from scales) are left alone:
  # they are small and refer to their own internal functions.
  in_frame <- function(e) {
    if (!is.environment(e) || identical(e, slim) || !is.null(attr(e, "class")) ||
        nzchar(environmentName(e))) {
      return(FALSE)
    }
    for (k in seq_len(1000)) {
      e <- parent.env(e)
      if (identical(e, emptyenv())) {
        return(FALSE)
      }
      nm <- environmentName(e)
      if (nzchar(nm)) {
        return(nm == home || startsWith(nm, "dartR"))
      }
    }
    FALSE
  }
  todo <- list() # names to copy, with the environment to look them up in

  queue <- function(names, env) {
    for (n in unique(names[nzchar(names)])) {
      todo[[length(todo) + 1]] <<- list(n = n, env = env)
    }
  }
  # a function defined inside the frame: refer to slim instead
  fix_fun <- function(f) {
    if (is.function(f) && !is.primitive(f) && in_frame(environment(f))) {
      queue(c(all.names(body(f)), unlist(lapply(formals(f), all.names))),
            environment(f))
      environment(f) <- slim
    }
    f
  }
  # an aesthetic mapping or facet variable (a quosure: a one-sided formula
  # with its environment as attribute)
  fix_quo <- function(q) {
    e <- attr(q, ".Environment")
    if (inherits(q, "formula") && in_frame(e)) {
      u <- unclass(q)
      queue(all.names(u[[length(u)]]), e)
      attr(q, ".Environment") <- slim
    }
    q
  }
  fix_list <- function(l) {
    if (is.list(l) && length(l)) {
      for (i in seq_along(l)) {
        if (inherits(l[[i]], "formula")) {
          l[[i]] <- fix_quo(l[[i]])
        } else if (is.function(l[[i]])) {
          l[[i]] <- fix_fun(l[[i]])
        }
      }
    }
    l
  }
  # is p a ggproto class defined in a package (e.g. ggplot2's
  # ScaleContinuousPosition), rather than an object made while plotting?
  is_class <- function(p) {
    cls <- class(p)[1]
    for (ns in loadedNamespaces()) {
      nse <- asNamespace(ns)
      if (exists(cls, envir = nse, inherits = FALSE) &&
          identical(get(cls, envir = nse, inherits = FALSE), p)) {
        return(TRUE)
      }
    }
    FALSE
  }
  # Flat copy of a ggproto object (an environment). ggproto objects inherit
  # fields through super(), a closure in the frame of the ggproto() call
  # whose unevaluated arguments lead back to the frame of the function that
  # made the plot; objects made while plotting (e.g. the copies ggplot2 makes
  # of scales when a plot is changed) inherit from each other this way. The
  # copy takes the fields of the object and of all such ancestors, and its
  # super() returns the first ancestor that is a class defined in a package,
  # so fetching a field gives the same value as before.
  clone <- function(o) {
    if (!is.environment(o)) {
      return(o)
    }
    fields <- list()
    e <- o
    parent_class <- NULL
    for (k in seq_len(100)) {
      f <- as.list.environment(e, all.names = TRUE)
      fields <- c(fields, f[setdiff(names(f), names(fields))])
      sup <- get0("super", envir = e, inherits = FALSE)
      if (!is.function(sup)) {
        break
      }
      p <- sup()
      if (!is.environment(p) || is_class(p)) {
        parent_class <- p
        break
      }
      e <- p
    }
    fields$super <- NULL
    o2 <- list2env(fields, parent = parent.env(o))
    attributes(o2) <- attributes(o)
    if (!is.null(parent_class)) {
      se <- new.env(parent = baseenv())
      assign("sv", parent_class, envir = se)
      assign("super", eval(quote(function() sv), se), envir = o2)
    }
    o2
  }

  slim_gg <- function(g) {
    if (inherits(g, "patchwork")) {
      for (i in seq_along(g)) {
        g[[i]] <- slim_gg(g[[i]])
      }
      return(g)
    }
    if (!inherits(g, "ggplot")) {
      return(g)
    }
    g$plot_env <- slim
    g$mapping <- fix_list(g$mapping)
    if (is.function(g$data)) {
      g$data <- fix_fun(g$data)
    }
    g$layers <- lapply(g$layers, function(l) {
      l <- clone(l)
      l$mapping <- fix_list(l$mapping)
      # set by ggplot2 when the plot has been printed (plot.display = TRUE)
      if (!is.null(l$computed_mapping)) {
        l$computed_mapping <- fix_list(l$computed_mapping)
      }
      if (is.function(l$data)) {
        l$data <- fix_fun(l$data)
      }
      l$aes_params <- fix_list(l$aes_params)
      l$geom_params <- fix_list(l$geom_params)
      l$stat_params <- fix_list(l$stat_params)
      l
    })
    if (is.environment(g$facet) && is.list(g$facet$params)) {
      fa <- clone(g$facet)
      pa <- fa$params
      for (nm in intersect(c("facets", "rows", "cols"), names(pa))) {
        pa[[nm]] <- fix_list(pa[[nm]])
      }
      fa$params <- pa
      g$facet <- fa
    }
    if (is.environment(g$scales) && is.list(g$scales$scales)) {
      sl <- clone(g$scales)
      sl$scales <- lapply(sl$scales, function(s) {
        if (!is.environment(s)) {
          return(s)
        }
        # always a copy: its super() may lead back to the frame (see clone())
        s <- clone(s)
        for (nm in ls(s, all.names = TRUE)) {
          f <- get(nm, envir = s)
          if (is.function(f)) {
            assign(nm, fix_fun(f), envir = s)
          }
        }
        s
      })
      g$scales <- sl
    }
    # the remaining ggproto parts: copies with a new super() (see clone())
    for (nm in c("layout", "coordinates", "guides")) {
      if (is.environment(g[[nm]])) {
        g[[nm]] <- clone(g[[nm]])
      }
    }
    g
  }

  out <- slim_gg(x)

  # copy the objects the plot refers to by name into slim; functions among
  # them that were defined in the frame are redirected to slim as well
  i <- 1
  while (i <= length(todo)) {
    n <- todo[[i]]$n
    e <- todo[[i]]$env
    i <- i + 1
    if (exists(n, envir = slim, inherits = FALSE)) {
      next
    }
    # look in the frame and the unnamed environments above it, up to anchor
    while (is.environment(e) && !identical(e, anchor) &&
           !identical(e, emptyenv()) && !identical(e, slim) &&
           !nzchar(environmentName(e))) {
      if (exists(n, envir = e, inherits = FALSE)) {
        assign(n, fix_fun(get(n, envir = e, inherits = FALSE)), envir = slim)
        break
      }
      e <- parent.env(e)
    }
  }
  out
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
