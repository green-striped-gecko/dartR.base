#' @name gl.plot.code
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
#' The code is printed at verbose >= 1 and returned invisibly. It is laid
#' out like a table in lines of at most 80 characters: code on the left, what
#' it does as a comment on the right, grouped in sections (titles and axes,
#' groups, fixed colours and sizes, single points, look, save). Each command
#' is runnable R code, so single commands can be copied. In the code,
#' the plot is called p (or by the name it was passed as). ggplot objects are
#' changed by adding layers, scales and themes with +, e.g.
#' p + labs(x = "Call rate"). In principle every component of a ggplot
#' (layers, scales, axes, labels, legend, theme) can be overwritten this way:
#' a component added with + replaces the one of the same kind already in the
#' plot. The code shown covers the most common changes; any other ggplot2
#' function can be added in the same way. A plot of several panels (a patchwork object)
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
#' gl.plot.code(p)
#'
#' @seealso \code{\link{gl.last.plot}}, \code{\link{theme_dartR}}
#' @export
#' @return The suggested code as a character vector, invisibly.

gl.plot.code <- function(p = NULL,
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
    code <- c(code, utils.plot.code.label("p <- gl.last.plot()",
                                          "last dartR plot"), "")
  }

  # FUNCTION SPECIFIC ERROR CHECKING
  # a recordedplot (base R graphics, see gl.last.plot()) is a list too, so it
  # must be recognised before the list branch
  is_plot <- function(z) {
    inherits(z, "ggplot") || inherits(z, "patchwork") || inherits(z, "recordedplot")
  }
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
      code <- c(code, utils.plot.code.one(p[[nm]], paste0(pname, "$", nm)), "")
    }
  } else if (is_plot(p)) {
    code <- c(code, utils.plot.code.one(p, pname), "")
  } else {
    stop(error("Fatal Error: p must be a ggplot or patchwork object, or a list",
               "of them as returned by gl.last.plot().\n"))
  }

  # DO THE JOB
  # lines of at most 80 characters: code on the left, what it does as a
  # comment on the right; every line stays runnable R code
  code <- c("# Every component of a ggplot can be overwritten by adding it with +;",
            "# the new one replaces the old one. Below are the most common changes.",
            "", code)
  if (verbose >= 1) {
    cat(code, sep = "\n")
  }
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  invisible(code)
}

# ---- layout of the printed code (lines of at most 80 characters) ----------

# a section header: "#-- title ------..." (or "#== title ===..." for a plot)
utils.plot.code.section <- function(title, fill = "-") {
  s <- paste0("#", fill, fill, " ", title, " ")
  if (nchar(s) > 76) s <- paste0(substr(s, 1, 72), "... ")
  paste0(s, strrep(fill, 80 - nchar(s)))
}

# code line(s) with what they do as a comment in a column on the right;
# if there is no room, the comment goes on its own line above the code
utils.plot.code.label <- function(lines, label) {
  first <- lines[1]
  col <- min(61, 80 - 3 - nchar(label))
  if (nchar(first) > col) {
    return(c(paste0("# ", label), lines))
  }
  c(paste0(sprintf(paste0("%-", col, "s"), first), " # ", label), lines[-1])
}

# head(arg1, arg2, ...) on one line, or wrapped after commas with the
# arguments aligned to the opening bracket
utils.plot.code.call <- function(head, args, limit = 61) {
  one <- paste0(head, "(", paste(args, collapse = ", "), ")")
  if (nchar(one) <= limit || length(args) == 1) return(one)
  pad <- strrep(" ", nchar(head) + 1)
  out <- character(0)
  cur <- paste0(head, "(")
  for (k in seq_along(args)) {
    piece <- paste0(args[k], if (k < length(args)) "," else ")")
    if (k > 1 && nchar(cur) + 1 + nchar(piece) > limit) {
      out <- c(out, cur)
      cur <- paste0(pad, piece)
    } else {
      cur <- paste0(cur, if (k > 1) " ", piece)
    }
  }
  c(out, cur)
}

# Code for one ggplot, patchwork or recordedplot object called pname
utils.plot.code.one <- function(p, pname) {
  if (inherits(p, "recordedplot")) {
    # base R graphics: no components to change, only redraw and save
    return(c(utils.plot.code.section(paste0(pname, ": base R graphics plot (recorded)"), "="),
             "# Not a ggplot: it cannot be changed with +. Change it with the arguments",
             "# of the function that drew it, then redraw or save it:",
             "",
             utils.plot.code.label(paste0("print(", pname, ")"), "redraw"),
             utils.plot.code.label("png(\"myplot.png\", width = 7, height = 5, units = \"in\", res = 300)",
                                   "open a file"),
             utils.plot.code.label(paste0("print(", pname, ")"), "draw into it"),
             utils.plot.code.label("dev.off()", "close the file")))
  }
  if (inherits(p, "patchwork")) {
    n <- length(p)
    code <- c(utils.plot.code.section(paste0(pname, ": figure of ", n,
                                             " panels (patchwork)"), "="),
              paste0("# change one panel with ", pname, "[[i]], all panels with &"),
              "",
              utils.plot.code.section("whole figure"),
              utils.plot.code.label(paste0(pname, " & theme(text = element_text(size = 14))"),
                                    "text size (all)"),
              utils.plot.code.label(paste0(pname, " + plot_annotation(title = \"My title\")"),
                                    "figure title"),
              utils.plot.code.save(pname))
    for (i in seq_len(n)) {
      sub <- p[[i]]
      if (!inherits(sub, "ggplot")) next
      code <- c(code, "", utils.plot.code.gg(sub, paste0(pname, "[[", i, "]]"),
                                             assign = TRUE))
    }
    return(code)
  }
  c(utils.plot.code.gg(p, pname, assign = FALSE), "",
    utils.plot.code.section("save"), utils.plot.code.save(pname))
}

# Code to save a plot
utils.plot.code.save <- function(pname) {
  utils.plot.code.label(paste0("ggsave(\"myplot.png\", ", pname,
                               ", width = 7, height = 5, dpi = 300)"),
                        "save (inches)")
}

# Code for one ggplot. With assign = TRUE (a panel of a patchwork) each line
# assigns the result back to the panel.
utils.plot.code.gg <- function(g, gname, assign = FALSE) {
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
  # "p + f" or, for a panel, "p[[i]] <- p[[i]] + f"
  addh <- function(f) {
    if (assign) paste0(gname, " <- ", gname, " + ", f) else paste0(gname, " + ", f)
  }
  call <- utils.plot.code.call
  label <- utils.plot.code.label
  section <- utils.plot.code.section

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

  # header: the plot, its title and the columns of its data
  ug <- unique(tolower(geoms))
  kind <- if (length(ug) == 0) {
    "empty plot"
  } else {
    paste0(ug[1], if (!grepl("plot$", ug[1])) " plot",
           if (length(ug) > 1) paste0(" (also ", paste(ug[-1], collapse = ", "), ")"))
  }
  out <- section(paste0(gname, ": ", kind,
                        if (!is.null(lab("title"))) paste0(" \"", lab("title"), "\"")), "=")
  if (is.data.frame(g$data) && ncol(g$data) > 0) {
    out <- c(out, strwrap(paste0("data ", gname, "$data, columns: ",
                                 paste(names(g$data), collapse = ", ")),
                          width = 78, prefix = "# "))
  }

  # titles and axes
  out <- c(out, "", section("titles and axes"),
           label(call(addh("labs"),
                      c(paste0("x = ", q(if (is.null(lab("x"))) "x title" else lab("x"))),
                        paste0("y = ", q(if (is.null(lab("y"))) "y title" else lab("y"))))),
                 "axis titles"),
           label(call(addh("labs"),
                      paste0("title = ", q(if (is.null(lab("title"))) "My title" else lab("title")))),
                 "plot title"))

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
    out <- c(out, label(call(addh("coord_cartesian"), lims), "axis ranges"))
  }

  # colours and symbols linked to a variable
  scale_lines <- function(aesthetic, fun_name, what) {
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
      pairs <- paste0(q(br[shown]), " = ", vtxt)
      lines <- call(paste0(addh(paste0("scale_", fun_name, "_manual")), "(values = c"),
                    pairs)
      lines[length(lines)] <- paste0(lines[length(lines)], ")")
      more <- if (length(br) > 10) {
        paste0("# (first 10 of ", length(br),
               " groups shown; add the others in the same way)")
      }
      c(label(lines, what), more)
    } else if (inherits(sc, "ScaleContinuous") && aesthetic %in% c("colour", "fill")) {
      label(call(addh(paste0("scale_", fun_name, "_gradient")),
                 c("low = \"yellow\"", "high = \"darkred\"")),
            paste(what, "(gradient)"))
    } else {
      NULL
    }
  }
  grp <- intersect(c("colour", "fill", "shape"), names(mapping))
  if (length(grp)) {
    vars <- unique(vapply(mapping[grp], aes_text, character(1)))
    gl_lines <- c(if ("colour" %in% grp) scale_lines("colour", "colour", "colours"),
                  if ("fill" %in% grp) scale_lines("fill", "fill", "fill colours"),
                  if ("shape" %in% grp) scale_lines("shape", "shape", "symbols"))
    if (length(gl_lines)) {
      out <- c(out, "", section(paste0("groups (", paste(vars, collapse = ", "), ")")),
               gl_lines)
    }
  }

  # constant colours, fills and sizes set in the layers
  fixed <- character(0)
  for (i in seq_along(g$layers)) {
    ap <- g$layers[[i]]$aes_params
    ap <- ap[intersect(names(ap), c("colour", "fill", "size", "linewidth", "alpha", "shape"))]
    for (a in names(ap)) {
      v <- ap[[a]]
      if (length(v) != 1) next
      vt <- if (is.character(v)) q(v) else as.character(v)
      fixed <- c(fixed, label(paste0(gname, "$layers[[", i, "]]$aes_params$", a, " <- ", vt),
                              paste0(a, ", layer ", i)))
    }
  }
  if (length(fixed)) {
    out <- c(out, "", section("fixed colours and sizes (change the plot in place)"), fixed)
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
      out <- c(out, "", section(paste0("single points (", idcol, ")")),
               label(paste0("one <- subset(", dname, ", ", idcol, " == ", q(val), ")"),
                     "select"),
               label(call(addh("geom_point"),
                          c("data = one", aes_xy, "colour = \"red\"", "size = 4",
                            "inherit.aes = FALSE")), "mark"),
               label(call(addh("geom_text"),
                          c("data = one",
                            paste0("aes(x = ", x, ", y = ", y, ", label = ", idcol, ")"),
                            "vjust = -1", "inherit.aes = FALSE")), "label"))
    } else {
      out <- c(out, "", section("single points (rows of the data)"),
               label(paste0("one <- ", dname, "[1, ]"), "select"),
               label(call(addh("geom_point"),
                          c("data = one", aes_xy, "colour = \"red\"", "size = 4",
                            "inherit.aes = FALSE")), "mark"),
               label(call(addh("geom_text"),
                          c("data = one", aes_xy, "label = \"my text\"",
                            "vjust = -1", "inherit.aes = FALSE")), "label"))
    }
  }

  # text and legend
  out <- c(out, "", section("look"),
           label(addh("theme(text = element_text(size = 14))"), "text size"),
           label(addh("theme(legend.position = \"bottom\")"), "legend"),
           label(addh("theme_dartR()"), "dartR theme"))
  out
}
