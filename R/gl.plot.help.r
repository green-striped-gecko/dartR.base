#' @name gl.plot.help
#' @title Shows the code to change a ggplot made by dartR
#' @family graphics
#'
#' @description
#' Looks at a ggplot (by default the plot of the last dartR function, see
#' \code{\link{gl.last.plot}}) and prints ready-to-use code to change it:
#' axis titles and ranges, the title, colours and symbols of the groups in
#' the plot, the colour of bars or lines, marking and labelling single data
#' points, text size, legend position and saving the plot. The code uses the
#' labels, groups and column names of the plot, so it can be copied, edited
#' and run.
#'
#' @param p A ggplot or patchwork object, or a named list of them as returned
#' by gl.last.plot(). If NULL, the plot(s) of the last dartR function
#' [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @details
#' The code is printed at verbose >= 1 and returned invisibly. In the code,
#' the plot is called p (or by the name it was passed as). ggplot objects are
#' changed by adding layers, scales and themes with +, e.g.
#' p + labs(x = "Call rate"). A plot of several panels (a patchwork object)
#' is changed per panel with p[[i]] <- p[[i]] + ..., or for all panels with
#' p & theme(...).
#'
#' Constant colours of bars or lines (not linked to a variable) are set in
#' the plot layer itself; the code shown changes them in place with
#' p$layers[[i]]$aes_params, which also changes any copy of p.
#'
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' pc <- gl.pcoa(possums.gl[1:50, 1:100], verbose = 0)
#' p <- gl.pcoa.plot(pc, possums.gl[1:50, 1:100], verbose = 0)
#' gl.plot.help(p)
#'
#' @seealso \code{\link{gl.last.plot}}, \code{\link{theme_dartR}}
#' @export
#' @return The suggested code as a character vector, invisibly.

gl.plot.help <- function(p = NULL,
                         verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  funname <- utils.flag.start(func = funname,
                              build = "v.2026.1",
                              verbose = verbose)

  # name of the plot in the code: the variable passed, otherwise p
  pname <- if (is.symbol(substitute(p))) deparse(substitute(p)) else "p"
  code <- character(0)
  if (is.null(p)) {
    p <- gl.last.plot(verbose = 0)
    pname <- "p"
    code <- c(code, "# the plot(s) of the last dartR function",
              "p <- gl.last.plot()")
  }

  # FUNCTION SPECIFIC ERROR CHECKING
  is_plot <- function(z) inherits(z, "ggplot") || inherits(z, "patchwork")
  if (is.list(p) && !is_plot(p)) {
    keep <- vapply(p, is_plot, logical(1))
    if (!any(keep)) {
      stop(error("Fatal Error: p contains no ggplot or patchwork object.\n"))
    }
    if (is.null(names(p)) || any(!nzchar(names(p)))) {
      names(p) <- paste0("plot", seq_along(p))
    }
    if (any(!keep) && verbose >= 2) {
      cat(report("  Not a plot, skipped:", paste0("$", names(p)[!keep], collapse = ", "), "\n"))
    }
    for (nm in names(p)[keep]) {
      code <- c(code, "", paste0("# ==== ", pname, "$", nm, " ===="),
                utils.plot.help.one(p[[nm]], paste0(pname, "$", nm)))
    }
  } else if (is_plot(p)) {
    code <- c(code, utils.plot.help.one(p, pname))
  } else {
    stop(error("Fatal Error: p must be a ggplot or patchwork object, or a list",
               "of them as returned by gl.last.plot().\n"))
  }

  # DO THE JOB
  if (verbose >= 1) {
    cat(code, sep = "\n")
    cat("\n")
  }
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  invisible(code)
}

# Code suggestions for one ggplot or patchwork object called pname
utils.plot.help.one <- function(p, pname) {
  if (inherits(p, "patchwork")) {
    n <- length(p)
    code <- c(paste0("# ", pname, " is a figure of ", n,
                     " panels (patchwork): change one panel with ",
                     pname, "[[i]], all panels with &"),
              paste0(pname, " & theme(text = element_text(size = 14))"),
              paste0(pname, " + plot_annotation(title = \"My title\")"))
    for (i in seq_len(n)) {
      sub <- p[[i]]
      if (!inherits(sub, "ggplot")) next
      code <- c(code, "", paste0("# ---- panel ", i, " ----"),
                utils.plot.help.gg(sub, paste0(pname, "[[", i, "]]"),
                                   assign = TRUE))
    }
    return(c(code, "", utils.plot.help.save(pname)))
  }
  c(utils.plot.help.gg(p, pname, assign = FALSE), "",
    utils.plot.help.save(pname))
}

# Code to save a plot
utils.plot.help.save <- function(pname) {
  c("# save the plot as an image (size in inches)",
    paste0("ggsave(\"myplot.png\", ", pname,
           ", width = 7, height = 5, dpi = 300)"))
}

# Code suggestions for one ggplot. With assign = TRUE (a panel of a
# patchwork) each line assigns the result back to the panel.
utils.plot.help.gg <- function(g, gname, assign = FALSE) {
  q <- function(s) paste0("\"", gsub("\"", "\\\\\"", s), "\"")
  # the expression of an aesthetic mapping: a quosure is a one-sided formula;
  # unclass() avoids rlang's deprecated [[ method for quosures
  aes_text <- function(m) {
    if (is.null(m)) return(NULL)
    e <- if (inherits(m, "formula")) {
      u <- unclass(m)
      u[[length(u)]]
    } else {
      m
    }
    paste(deparse(e), collapse = "")
  }
  add <- function(x) {
    if (assign) paste0(gname, " <- ", gname, " + ", x) else paste0(gname, " + ", x)
  }
  out <- character(0)

  labs_now <- tryCatch(
    if (exists("get_labs", envir = asNamespace("ggplot2"))) {
      ggplot2::get_labs(g)
    } else {
      g$labels
    }, error = function(e) list())
  lab <- function(a) {
    v <- labs_now[[a]]
    if (is.character(v) && length(v) == 1) v else NULL
  }

  # aesthetics: those of the plot, completed by those of its layers
  mapping <- as.list(g$mapping)
  geoms <- character(0)
  for (l in g$layers) {
    geoms <- c(geoms, sub("^Geom", "", class(l$geom)[1]))
    for (a in names(l$mapping)) {
      if (is.null(mapping[[a]])) mapping[[a]] <- l$mapping[[a]]
    }
  }
  x <- aes_text(mapping$x)
  y <- aes_text(mapping$y)

  out <- c(out, paste0("# ", if (length(geoms)) paste(unique(geoms), collapse = ", ") else "no layers",
                       if (!is.null(lab("title"))) paste0(": \"", lab("title"), "\"") else ""))

  # titles
  labs_args <- c(paste0("x = ", q(if (is.null(lab("x"))) "x title" else lab("x"))),
                 paste0("y = ", q(if (is.null(lab("y"))) "y title" else lab("y"))),
                 paste0("title = ", q(if (is.null(lab("title"))) "My title" else lab("title"))))
  out <- c(out, "# axis titles and plot title",
           add(paste0("labs(", paste(labs_args, collapse = ", "), ")")))

  # axis ranges of the drawn plot
  built <- tryCatch(suppressMessages(suppressWarnings(ggplot2::ggplot_build(g))),
                    error = function(e) NULL)
  rng <- function(axis) {
    pp <- tryCatch(built$layout$panel_params[[1]][[paste0(axis, ".range")]],
                   error = function(e) NULL)
    if (is.numeric(pp) && length(pp) == 2 && all(is.finite(pp))) {
      signif(pp, 3)
    } else {
      NULL
    }
  }
  xr <- rng("x")
  yr <- rng("y")
  if (!is.null(xr) || !is.null(yr)) {
    lims <- c(if (!is.null(xr)) paste0("xlim = c(", xr[1], ", ", xr[2], ")"),
              if (!is.null(yr)) paste0("ylim = c(", yr[1], ", ", yr[2], ")"))
    out <- c(out, "# axis ranges (zoom in or out; the values are the current ranges)",
             add(paste0("coord_cartesian(", paste(lims, collapse = ", "), ")")))
  }

  # colours and symbols linked to a variable
  scale_lines <- function(aesthetic, fun_name) {
    sc <- tryCatch(built$plot$scales$get_scales(aesthetic), error = function(e) NULL)
    if (is.null(sc)) return(NULL)
    if (inherits(sc, "ScaleDiscrete")) {
      br <- tryCatch(sc$get_breaks(), error = function(e) NULL)
      br <- br[!is.na(br)]
      if (length(br) == 0) return(NULL)
      vals <- tryCatch(sc$map(br), error = function(e) NULL)
      if (is.null(vals) || length(vals) != length(br)) return(NULL)
      shown <- seq_len(min(length(br), 10))
      vtxt <- if (is.numeric(vals)) as.character(vals[shown]) else q(vals[shown])
      body <- paste0(q(br[shown]), " = ", vtxt, collapse = ", ")
      more <- if (length(br) > 10) {
        paste0("# (the first 10 of ", length(br),
               " groups are shown; add the others in the same way)")
      }
      c(paste0("# ", aesthetic, " of each group of ", aes_text(mapping[[aesthetic]]),
               " (the values are the current ones)"),
        more,
        add(paste0("scale_", fun_name, "_manual(values = c(", body, "))")))
    } else if (inherits(sc, "ScaleContinuous") && aesthetic %in% c("colour", "fill")) {
      c(paste0("# ", aesthetic, " gradient of ", aes_text(mapping[[aesthetic]])),
        add(paste0("scale_", fun_name, "_gradient(low = \"yellow\", high = \"darkred\")")))
    } else {
      NULL
    }
  }
  if (!is.null(mapping$colour)) out <- c(out, scale_lines("colour", "colour"))
  if (!is.null(mapping$fill)) out <- c(out, scale_lines("fill", "fill"))
  if (!is.null(mapping$shape)) out <- c(out, scale_lines("shape", "shape"))

  # constant colours, fills and sizes set in the layers
  for (i in seq_along(g$layers)) {
    ap <- g$layers[[i]]$aes_params
    ap <- ap[intersect(names(ap), c("colour", "fill", "size", "linewidth", "alpha", "shape"))]
    if (length(ap) == 0) next
    out <- c(out, paste0("# fixed settings of layer ", i, " (", geoms[i],
                         "): changes the plot in place"))
    for (a in names(ap)) {
      v <- ap[[a]]
      if (length(v) != 1) next
      vt <- if (is.character(v)) q(v) else as.character(v)
      out <- c(out, paste0(gname, "$layers[[", i, "]]$aes_params$", a, " <- ", vt))
    }
  }

  # marking single points
  if (any(geoms %in% c("Point", "Text", "Label")) && !is.null(x) && !is.null(y) &&
      is.data.frame(g$data) && nrow(g$data) > 0) {
    d <- g$data
    used <- vapply(mapping, aes_text, character(1))
    chr <- names(d)[vapply(d, function(z) is.character(z) || is.factor(z), logical(1))]
    idcol <- c(intersect(c("ind.names", "ind", "id", "label", "labels", "locus", "loc.names"), chr),
               setdiff(chr, used), chr)[1]
    dname <- paste0(gname, "$data")
    aes_xy <- paste0("aes(x = ", x, ", y = ", y, ")")
    if (!is.na(idcol)) {
      val <- as.character(d[[idcol]][1])
      sel <- paste0("subset(", dname, ", ", idcol, " == ", q(val), ")")
      out <- c(out,
               paste0("# mark single points: here those with ", idcol, " == ", q(val),
                      " (columns of the data: ", paste(names(d), collapse = ", "), ")"),
               add(paste0("geom_point(data = ", sel, ", ", aes_xy,
                          ", colour = \"red\", size = 4, inherit.aes = FALSE)")),
               "# and label them",
               add(paste0("geom_text(data = ", sel, ", aes(x = ", x, ", y = ", y,
                          ", label = ", idcol, "), vjust = -1, inherit.aes = FALSE)")))
    } else {
      out <- c(out,
               paste0("# mark single points: select rows of the data (columns: ",
                      paste(names(d), collapse = ", "), ")"),
               add(paste0("geom_point(data = ", dname, "[1, ], ", aes_xy,
                          ", colour = \"red\", size = 4, inherit.aes = FALSE)")))
    }
  }

  # text and legend
  out <- c(out, "# text size, legend position, theme",
           add("theme(text = element_text(size = 14), legend.position = \"bottom\")"),
           add("theme_dartR()"))
  out
}
