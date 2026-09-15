#' @name gl.plot.heatmap
#' @title Represents a distance matrix as a heatmap
#' @family graphics

#' @description
#' Plots a heat map of the values in a distance or dissimilarity matrix,
#' with optional row and column dendrograms. When a genlight object is
#' supplied, individuals are coloured by population in side bars and a
#' legend. The heat map is drawn by \code{\link{utils.heatmap}}, a modified
#' copy of \code{gplots::heatmap.2}.
#'
#' @details
#' \code{D} can be a \code{dist} object (for example from
#' \code{\link{gl.dist.ind}} or \code{\link{gl.dist.pop}}), a square numeric
#' matrix with matching row and column names, or an object of class
#' \code{fd} from \code{\link{gl.fixed.diff}}, in which case the raw fixed
#' difference matrix (\code{$fd}) is plotted. A matrix is plotted as
#' supplied: its diagonal and both triangles are drawn, so a relatedness
#' matrix keeps its self-relatedness values and an asymmetric matrix keeps
#' both directions. If one triangle is entirely \code{NA}, it is filled from
#' the other. Set \code{diag.na = TRUE} to leave the diagonal blank.
#'
#' Population colours are drawn when every column of \code{D} names an
#' individual of \code{x}. When \code{D} is at the population level
#' (\code{gl.dist.pop}, \code{fd}), \code{x} is ignored. When some columns
#' of \code{D} are not individuals of \code{x}, no colours are drawn and a
#' warning is printed at \code{verbose >= 1}.
#'
#' \code{legendx} and \code{legendy} position the legend in the 0-1
#' coordinates of the plotting region; \code{legendx} also accepts a
#' keyword such as \code{"topleft"} (see \code{\link[graphics]{legend}}).
#'
#' @param D Distance matrix (class \code{dist}), square numeric matrix, or
#' object of class \code{fd} [required].
#' @param x Genlight object used to colour individuals by population
#' [default NULL].
#' @param palette.divergent A divergent palette function for the values
#' [default gl.colors("div")].
#' @param palette_discrete Colours for the populations of \code{x}: a
#' palette function, or a vector with one colour per population
#' [default NULL, in which case colours are chosen by gl.select.colors].
#' @param dendrogram Which dendrograms to draw: 'none', 'row', 'column' or
#' 'both' [default "column"].
#' @param plot.out Whether to draw the heat map [default TRUE].
#' @param legend.print Whether to draw the population legend (only when
#' \code{x} is supplied and matches \code{D}) [default TRUE].
#' @param legendx x position of the legend, 0-1, or a legend keyword
#' [default 0].
#' @param legendy y position of the legend, 0-1 [default 0.5].
#' @param label.size Size of the population labels in the legend
#' [default 0.75].
#' @param legend.title Legend title [default "Populations"].
#' @param diag.na If TRUE, the diagonal of the matrix is set to NA and drawn
#' in \code{na.color} [default FALSE].
#' @param margins Numeric vector of length 2: margins for column and row
#' names, respectively [default c(10, 10)].
#' @param na.color Colour for missing values (NA) [default "grey"].
#' @param revC Reverse column order [default FALSE].
#' @param symbreaks Symmetric colour breaks around zero [default FALSE].
#' @param trace Trace lines: "column", "row", "both" or "none"
#' [default "none"].
#' @param tracecol Trace line colour [default "cyan"].
#' @param cexRow Row label scale [default NULL, chosen from the number of
#' rows].
#' @param cexCol Column label scale [default NULL, chosen from the number
#' of columns].
#' @param srtRow Row label rotation angle [default NULL].
#' @param srtCol Column label rotation angle [default 90].
#' @param offsetRow Row label offset [default 0.5].
#' @param offsetCol Column label offset [default 0.5].
#' @param key Draw the colour key [default TRUE].
#' @param keysize Size of the colour key [default 1.5].
#' @param density.info Density plot on the colour key: "histogram",
#' "density" or "none" [default "none"].
#' @param denscol Density plot colour [default "cyan"].
#' @param symkey Symmetric colour key around zero [default FALSE].
#' @param densadj Density plot bandwidth adjustment [default 0.25].
#' @param key.title Title of the colour key [default NULL].
#' @param key.xlab X label of the colour key [default NULL].
#' @param key.ylab Y label of the colour key [default NULL].
#' @param main Main plot title [default NULL].
#' @param xlab X-axis label [default NULL].
#' @param ylab Y-axis label [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @param ... Further arguments passed to \code{\link{utils.heatmap}}.
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' if (requireNamespace("dendextend", quietly = TRUE)) {
#'   gl <- testset.gl[1:10, ]
#'   D <- dist(as.matrix(gl))
#'   gl.plot.heatmap(D, x = gl)
#'   D2 <- gl.dist.pop(possums.gl)
#'   gl.plot.heatmap(D2)
#' }
#' \donttest{
#' D3 <- gl.fixed.diff(testset.gl)
#' gl.plot.heatmap(D3)
#' }
#' @importFrom graphics legend par
#' @export
#' @return Invisibly, the list returned by \code{\link{utils.heatmap}}
#' (row and column orders, the plotted matrix, dendrograms, breaks and
#' colours), or NULL when \code{plot.out = FALSE}.

gl.plot.heatmap <- function(D,
                            x = NULL,
                            palette.divergent = gl.colors("div", verbose = 0),
                            palette_discrete = NULL,
                            dendrogram = "column",
                            plot.out = TRUE,
                            legend.print = TRUE,
                            legendx = 0,
                            legendy = 0.5,
                            label.size = 0.75,
                            legend.title = "Populations",
                            diag.na = FALSE,
                            margins = c(10, 10),
                            na.color = "grey",
                            revC = FALSE,
                            symbreaks = FALSE,
                            trace = "none",
                            tracecol = "cyan",
                            cexRow = NULL,
                            cexCol = NULL,
                            srtRow = NULL,
                            srtCol = 90,
                            offsetRow = 0.5,
                            offsetCol = 0.5,
                            key = TRUE,
                            keysize = 1.5,
                            density.info = "none",
                            denscol = "cyan",
                            symkey = FALSE,
                            densadj = 0.25,
                            key.title = NULL,
                            key.xlab = NULL,
                            key.ylab = NULL,
                            main = NULL,
                            xlab = NULL,
                            ylab = NULL,
                            verbose = NULL,
                            ...) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <-
    utils.check.datatype(D,
                         accept = c("dist", "fd", "matrix"),
                         verbose = verbose)
  
  # CHECK IF PACKAGES ARE INSTALLED
  pkg <- "dendextend"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
  }
  
  # FUNCTION SPECIFIC ERROR CHECKING
  if (!is.null(x) && !is.null(palette_discrete) &&
      !is.function(palette_discrete) &&
      length(palette_discrete) != nPop(x)) {
    stop(error(
      "palette_discrete must be a palette function or a vector of",
      nPop(x),
      "colours, one per population of x; got",
      length(palette_discrete),
      "\n"
    ))
  }
  
  # DO THE JOB
  
  # The matrix to draw
  if (datatype == "dist") {
    m <- as.matrix(D)
  } else if (datatype == "matrix") {
    m <- D
    if (nrow(m) != ncol(m)) {
      stop(error(
        "D must be a square matrix; got",
        nrow(m), "rows and", ncol(m), "columns\n"
      ))
    }
    if (!is.null(rownames(m)) && !is.null(colnames(m)) &&
        !identical(rownames(m), colnames(m))) {
      stop(error("Row and column names of D must be identical\n"))
    }
    # A matrix with one triangle left empty is mirrored from the other
    up <- upper.tri(m)
    lo <- lower.tri(m)
    if (all(is.na(m[up])) && !all(is.na(m[lo]))) {
      m[up] <- t(m)[up]
    } else if (all(is.na(m[lo])) && !all(is.na(m[up]))) {
      m[lo] <- t(m)[lo]
    }
  } else {
    m <- as.matrix(D$fd)
  }
  
  if (diag.na) {
    diag(m) <- NA
  }
  
  # Population colours for the side bars, labels and legend
  colors_pops <- NULL
  legend_text <- NULL
  legend_color <- NULL
  
  if (is.null(x)) {
    legend.print <- FALSE
  } else {
    if (is.null(palette_discrete)) {
      pop_colors <- gl.select.colors(x, verbose = 0)
    } else if (is.function(palette_discrete)) {
      pop_colors <- palette_discrete(nPop(x))
    } else {
      pop_colors <- palette_discrete
    }
    names(pop_colors) <- as.character(popNames(x))
    
    idx <- match(colnames(m), indNames(x))
    if (length(idx) == ncol(m) && !anyNA(idx)) {
      # every column of D is an individual of x
      ind_pops <- as.character(pop(x))[idx]
      colors_pops <- unname(pop_colors[ind_pops])
      legend_text <- popNames(x)[popNames(x) %in% ind_pops]
      legend_color <- unname(pop_colors[legend_text])
    } else if (datatype == "fd" ||
               (!is.null(colnames(m)) && all(colnames(m) %in% popNames(x)))) {
      # D is at the population level; nothing to colour
      if (verbose >= 2) {
        cat(report("  D is a population-level matrix; x is ignored\n"))
      }
      legend.print <- FALSE
    } else {
      n_missing <- if (is.null(colnames(m))) ncol(m) else sum(is.na(idx))
      if (verbose >= 1) {
        cat(warn(
          "  Warning:", n_missing, "of", ncol(m),
          "columns of D are not individuals of x;",
          "population colours and legend not drawn\n"
        ))
      }
      legend.print <- FALSE
    }
  }
  
  p3 <- NULL
  if (plot.out) {
    p3 <-
      utils.heatmap(
        m,
        col = palette.divergent(255),
        dendrogram = dendrogram,
        margins = margins,
        na.color = na.color,
        ColSideColors = colors_pops,
        RowSideColors = colors_pops,
        colRow = colors_pops,
        colCol = colors_pops,
        revC = revC,
        symbreaks = symbreaks,
        trace = trace,
        tracecol = tracecol,
        cexRow = cexRow,
        cexCol = cexCol,
        srtRow = srtRow,
        srtCol = srtCol,
        offsetRow = offsetRow,
        offsetCol = offsetCol,
        key = key,
        keysize = keysize,
        density.info = density.info,
        denscol = denscol,
        symkey = symkey,
        densadj = densadj,
        key.title = key.title,
        key.xlab = key.xlab,
        key.ylab = key.ylab,
        main = main,
        xlab = xlab,
        ylab = ylab,
        ...
      )
    if (legend.print) {
      op <- par(mar = c(1, 1, 1, 1))
      on.exit(par(op), add = TRUE)
      legend(
        legendx,
        legendy,
        legend = legend_text,
        fill = legend_color,
        cex = label.size,
        title = legend.title
      )
    }
  }
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  invisible(p3)
}
