#' @name utils.heatmap
#' @title Enhanced Heat Map
#'
#' @description
#' A heat map is a false‐color image with optional row/column dendrograms
#' and extensive customization of scaling, coloring, labeling, and layout.
#'
#' @param x Numeric matrix of the values to be plotted.
#' @param Rowv Logical, dendrogram, or vector. Controls row‐dendrogram creation and/or reordering.
#' @param Colv Logical, dendrogram, or vector. Controls column‐dendrogram creation and/or reordering; may be "Rowv" for square matrices.
#' @param distfun Function to compute distances; defaults to \code{\link{dist}}.
#' @param hclustfun Function to cluster; defaults to \code{\link{hclust}}.
#' @param dendrogram Which dendrogram(s) to draw: "both", "row", "column", or "none".
#' @param reorderfun Function to reorder dendrogram given weights.
#' @param symm Logical; treat \code{x} symmetrically (only for square matrices).
#' @param scale Character: "none", "row", or "column" scaling of \code{x}.
#' @param na.rm Logical; remove \code{NA}s in scaling and computations.
#' @param revC Logical; reverse column order for plotting.
#' @param add.expr Expression to evaluate after the heat map image is drawn.
#' @param breaks Numeric vector or integer: breakpoints for color bins.
#' @param symbreaks Logical; force symmetric breaks around zero.
#' @param col Color palette or function; defaults to \code{heat.colors}.
#' @param colsep,rowsep Integer vectors; where to put separator lines between blocks.
#' @param sepcolor Color for separators.
#' @param sepwidth Numeric vector of length 2; widths of separators relative to cell size.
#' @param cellnote Optional character matrix for annotating each cell.
#' @param notecex Numeric; scaling for \code{cellnote}.
#' @param notecol Color for \code{cellnote}; default "cyan".
#' @param na.color Color for \code{NA} cells; defaults to plotting background.
#' @param trace Character: "column", "row", "both", or "none" for trace lines.
#' @param tracecol Color for trace lines.
#' @param hline,vline Numeric vectors; values at which to draw horizontal/vertical lines.
#' @param linecol Color for \code{hline}/\code{vline}.
#' @param margins Numeric vector of length 2; margins for column and row labels.
#' @param ColSideColors,RowSideColors Optional color vectors for side bars [default NULL].
#' @param cexRow,cexCol Numeric; character expansion for row/column labels.
#' @param labRow,labCol Character vectors of labels; defaults to row/column names.
#' @param srtRow,srtCol Rotation angles for row/column labels.
#' @param adjRow,adjCol Justification for row/column labels.
#' @param offsetRow,offsetCol Numeric offsets (in character widths) for labels.
#' @param colRow,colCol Color(s) for row/column label text.
#' @param key Logical; draw a color key.
#' @param keysize Numeric; size of the color key.
#' @param density.info Character: "histogram", "density", or "none" on the key.
#' @param denscol Color for density plot.
#' @param symkey Logical; symmetric key around zero.
#' @param densadj Numeric; adjustment for density bandwidth.
#' @param key.title,key.xlab,key.ylab Main title and axis labels for key.
#' @param key.xtickfun,key.ytickfun Functions for tick placement on key axes.
#' @param key.par List of graphical parameters for the key.
#' @param main,xlab,ylab Main title and axis labels for the heat map.
#' @param lmat Matrix specifying layout positions.
#' @param lhei Numeric vector specifying row heights in layout.
#' @param lwid Numeric vector specifying column widths in layout.
#' @param extrafun Function to call after plotting (e.g. add scatterplot).
#' @param ... Additional arguments passed to \code{\link{image}}.
#'
#' @details
#' If \code{Rowv} or \code{Colv} are dendrograms, they are honored
#' without reordering; otherwise they are computed via
#' \code{as.dendrogram(hclustfun(distfun(...)))} and optionally reordered.
#'
#' Scaling centers and scales rows or columns when requested.  The
#' default color palette may be enhanced via packages such as
#' \pkg{RColorBrewer}.  Layout may be customized with
#' \code{lmat}, \code{lwid}, and \code{lhei} as in \code{\link[graphics]{layout}}.
#'
#' @note
#' The original rows and columns are reordered to match the dendrograms.
#' Because \code{heatmap.2} uses \code{\link{layout}}, it cannot be used
#' inside another layout (e.g. \code{par(mfrow=...)}).
#'
#' @return Invisibly, a list with components:
#' \describe{
#'   \item{rowInd}{Row permutation vector.}
#'   \item{colInd}{Column permutation vector.}
#'   \item{call}{Matched function call.}
#'   \item{rowMeans,rowSDs}{Row means and SDs (if \code{scale="row"}).}
#'   \item{colMeans,colSDs}{Column means and SDs (if \code{scale="column"}).}
#'   \item{carpet}{Reordered/scaled matrix used for the image.}
#'   \item{rowDendrogram,colDendrogram}{Dendrogram objects.}
#'   \item{breaks}{Breakpoints used for colors.}
#'   \item{col}{Colors used.}
#'   \item{vline,hline}{Center‐line values for trace (if used).}
#'   \item{colorTable}{Data frame of color bins and ranges.}
#'   \item{layout}{List with \code{lmat}, \code{lhei}, \code{lwid}.}
#' }
#'
#' @author
#' Andy Liaw (original); R. Gentleman, M. Maechler, W. Huber, G. Warnes (revisions)
#'
#' @importFrom graphics layout strheight strwidth mtext rect abline plot.new title
#'
#' @examples
#' data(mtcars)
#' x <- as.matrix(mtcars)
#' utils.heatmap(x)
#' utils.heatmap(x, dendrogram = "none")
#' utils.heatmap(x, scale = "row", col = cm.colors(255))
#' @export

utils.heatmap <- function(x,
                      Rowv = TRUE,
                      Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both", "row", "column", "none"),
                      reorderfun = function(d, w) reorder(d, w),
                      symm = FALSE,
                      scale = c("none", "row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv, "Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = any(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("column", "row", "both", "none"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5, 5),
                      ColSideColors = NULL,
                      RowSideColors = NULL,
                      cexRow = 0.2 + 1 / log10(nr),
                      cexCol = 0.2 + 1 / log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      srtRow = NULL,
                      srtCol = NULL,
                      adjRow = c(0, NA),
                      adjCol = c(NA, 0),
                      offsetRow = 0.5,
                      offsetCol = 0.5,
                      colRow = NULL,
                      colCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("histogram", "density", "none"),
                      denscol = tracecol,
                      symkey = any(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      key.title = NULL,
                      key.xlab = NULL,
                      key.ylab = NULL,
                      key.xtickfun = NULL,
                      key.ytickfun = NULL,
                      key.par = list(),
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      extrafun = NULL,
                      ...) {
  
  # Helper function to scale values between 0 and 1
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low) / (high - low)
    x
  }
  
  margin_dendrogram <- (max(nchar(colnames(x)))/2) + 0.5
  
  # Initialize return value list
  retval <- list()
  
  # Match arguments
  scale <- if (symm && missing(scale)) "none" else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  
  # Check color function
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  
  # Validate breaks
  if (!missing(breaks) && any(duplicated(breaks)))
    stop("breaks may not contain duplicate values")
  
  if (!missing(breaks) && (scale != "none"))
    warning(
      "Using scale=\"row\" or scale=\"column\" when breaks are ",
      "specified can produce unpredictable results. ",
      "Please consider using only one or the other."
    )
  
  # Validate row and column specifications
  if (is.null(Rowv) || any(is.na(Rowv)))
    Rowv <- FALSE
  
  if (is.null(Colv) || any(is.na(Colv)))
    Colv <- FALSE
  else if (all(Colv == "Rowv"))
    Colv <- Rowv
  
  # Check input matrix
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("x must be a numeric matrix")
  
  nr <- di[1]
  nc <- di[2]
  
  if (nr <= 1 || nc <= 1)
    stop("x must have at least 2 rows and 2 columns")
  
  if (!is.numeric(margins) || length(margins) != 2)
    stop("margins must be a numeric vector of length 2")
  
  # Initialize cell notes if missing
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  
  # Check for dendrogram discrepancies
  if (!inherits(Rowv, "dendrogram")) {
    if (((is.logical(Rowv) && !isTRUE(Rowv)) || (is.null(Rowv))) &&
        (dendrogram %in% c("both", "row"))) {
      warning(
        "Discrepancy: Rowv is FALSE, while dendrogram is ",
        dendrogram,
        ". Omitting row dendrogram."
      )
      if (dendrogram == "both")
        dendrogram <- "column"
      else
        dendrogram <- "none"
    }
  }
  
  if (!inherits(Colv, "dendrogram")) {
    if (((is.logical(Colv) && !isTRUE(Colv)) || (is.null(Colv))) &&
        (dendrogram %in% c("both", "column"))) {
      warning(
        "Discrepancy: Colv is FALSE, while dendrogram is ",
        dendrogram,
        ". Omitting column dendrogram."
      )
      if (dendrogram == "both")
        dendrogram <- "row"
      else
        dendrogram <- "none"
    }
  }
  
  # Handle row clustering and ordering
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
    if (length(rowInd) > nr || any(rowInd < 1 | rowInd > nr))
      stop("Rowv dendrogram doesn't match size of x")
    if (length(rowInd) < nr)
      nr <- length(rowInd)
  } else if (is.integer(Rowv)) {
    distr <- distfun(x)
    hcr <- hclustfun(distr)
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  } else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    distr <- distfun(x)
    hcr <- hclustfun(distr)
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  } else if (!isTRUE(Rowv)) {
    rowInd <- nr:1
    ddr <- as.dendrogram(hclust(dist(diag(nr))))
  } else {
    rowInd <- nr:1
    ddr <- as.dendrogram(Rowv)
  }
  
  # Handle column clustering and ordering
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
    if (length(colInd) > nc || any(colInd < 1 | colInd > nc))
      stop("Colv dendrogram doesn't match size of x")
    if (length(colInd) < nc)
      nc <- length(colInd)
  } else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    } else
      colInd <- rowInd
  } else if (is.integer(Colv)) {
    distc <- distfun(if (symm) x else t(x))
    hcc <- hclustfun(distc)
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  } else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    distc <- distfun(if (symm) x else t(x))
    hcc <- hclustfun(distc)
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  } else if (!isTRUE(Colv)) {
    colInd <- 1:nc
    ddc <- as.dendrogram(hclust(dist(diag(nc))))
  } else {
    colInd <- 1:nc
    ddc <- as.dendrogram(Colv)
  }
  
  # Store permutation vectors
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  
  # Reorder data
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  
  # Prepare labels
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else
    rownames(x)
  else
    labRow <- labRow[rowInd]
  
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else
    colnames(x)
  else
    labCol <- labCol[colInd]
  
  if (!is.null(colRow))
    colRow <- colRow[rowInd]
  
  if (!is.null(colCol))
    colCol <- colCol[colInd]
  
  # Scale data if requested
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  } else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  
  # Set up color breaks
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else
      breaks <- length(col) + 1
  }
  
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  
  if (is(col, "function"))
    col <- col(ncol)
  
  # Clip values to breaks range
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  
  # Set up layout dimensions
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  
  # Set up layout matrix
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!is.null(ColSideColors)) {
      if (!is.character(ColSideColors) || length(ColSideColors) != nc)
        stop("ColSideColors must be a character vector of length ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      lhei <- c(lhei[1], 0.2, lhei[2])
    }
    
    if (!is.null(RowSideColors)) {
      if (!is.character(RowSideColors) || length(RowSideColors) != nr)
        stop("RowSideColors must be a character vector of length nrow(x)")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], 0.2, lwid[2])
    }
    
    lmat[is.na(lmat)] <- 0
  }
  
  # Check layout dimensions
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  
  # Save graphics parameters and set up layout
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  # Plot index counter
  plot.index <- 1
  
  # Draw row side colors if requested
  if (!is.null(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0))
    image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    plot.index <- plot.index + 1
  }
  
  # Draw column side colors if requested
  if (!is.null(ColSideColors)) {
    par(mar = c(0, 0, 0, margins[2]))
    image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    plot.index <- plot.index + 1
  }
  
  # Main heatmap plot
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  } else {
    iy <- 1:nr
  }
  
  # Draw the heatmap
  image(
    1:nc,
    1:nr,
    x,
    xlim = 0.5 + c(0, nc),
    ylim = 0.5 + c(0, nr),
    axes = FALSE,
    xlab = "",
    ylab = "",
    col = col,
    breaks = breaks,
    ...
  )
  
  retval$carpet <- x
  
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  
  retval$breaks <- breaks
  retval$col <- col
  
  # Handle NA values
  
  invalid2 <-  function (x) 
  {
    if (missing(x) || is.null(x) || length(x) == 0 || is(x, "try-error")) {
      return(TRUE)
    }
    if (is.list(x)) {
      return(all(sapply(x, invalid2)))
    }
    else if (is.vector(x)) {
      return(all(is.na(x)))
    }
    else {
      return(FALSE)
    }
  }
  
  
  if (!invalid2(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(
      1:nc,
      1:nr,
      mmat,
      axes = FALSE,
      xlab = "",
      ylab = "",
      col = na.color,
      add = TRUE
    )
  }
  
  # Add column labels
  if (is.null(srtCol) && is.null(colCol)) {
    axis(
      1,
      1:nc,
      labels = labCol,
      las = 2,
      line = -0.5 + offsetCol,
      tick = 0,
      cex.axis = cexCol,
      hadj = adjCol[1],
      padj = adjCol[2]
    )
  } else {
    if (is.null(srtCol) || is.numeric(srtCol)) {
      if (missing(adjCol) || is.null(adjCol))
        adjCol = c(1, NA)
      
      if (is.null(srtCol))
        srtCol <- 90
      
      xpd.orig <- par("xpd")
      par(xpd = NA)
      
      xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, tick = 0)
      
      text(
        x = xpos,
        y = par("usr")[3] - (1 + offsetCol) * strheight("M"),
        labels = labCol,
        adj = adjCol,
        cex = cexCol,
        srt = srtCol,
        col = colCol
      )
      
      par(xpd = xpd.orig)
    } else {
      warning("Invalid value for srtCol ignored.")
    }
  }
  
  # Add row labels
  if (is.null(srtRow) && is.null(colRow)) {
    axis(
      4,
      iy,
      labels = labRow,
      las = 2,
      line = -0.5 + offsetRow,
      tick = 0,
      cex.axis = cexRow,
      hadj = adjRow[1],
      padj = adjRow[2]
    )
  } else {
    if (is.null(srtRow) || is.numeric(srtRow)) {
      xpd.orig <- par("xpd")
      par(xpd = NA)
      
      ypos <- axis(4, iy, labels = rep("", nr), las = 2, line = -0.5, tick = 0)
      
      text(
        x = par("usr")[2] + (1 + offsetRow) * strwidth("M"),
        y = ypos,
        labels = labRow,
        adj = adjRow,
        cex = cexRow,
        srt = srtRow,
        col = colRow
      )
      
      par(xpd = xpd.orig)
    } else {
      warning("Invalid value for srtRow ignored.")
    }
  }
  
  # Add axis labels
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  
  # Add custom expressions
  if (!missing(add.expr))
    eval(substitute(add.expr))
  
  # Add column separators
  if (!missing(colsep)) {
    for (csep in colsep) {
      rect(
        xleft = csep + 0.5,
        ybottom = 0,
        xright = csep + 0.5 + sepwidth[1],
        ytop = ncol(x) + 1,
        lty = 1,
        lwd = 1,
        col = sepcolor,
        border = sepcolor
      )
    }
  }
  
  # Add row separators
  if (!missing(rowsep)) {
    for (rsep in rowsep) {
      rect(
        xleft = 0,
        ybottom = (ncol(x) + 1 - rsep) - 0.5,
        xright = nrow(x) + 1,
        ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2],
        lty = 1,
        lwd = 1,
        col = sepcolor,
        border = sepcolor
      )
    }
  }
  
  # Scale for trace lines
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  
  # Add column trace lines
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    
    for (i in 1:length(colInd)) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol, lty = 2)
      }
      
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  
  # Add row trace lines
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    
    for (i in 1:length(rowInd)) {
      if (!is.null(hline)) {
        abline(h = i - 0.5 + hline.vals, col = linecol, lty = 2)
      }
      
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  
  # Add cell annotations
  if (!missing(cellnote)) {
    text(
      x = c(row(cellnote)),
      y = c(col(cellnote)),
      labels = c(cellnote),
      col = notecol,
      cex = notecex
    )
  }
  
  plot.index <- plot.index + 1
  
  # Draw row dendrogram
  # par(mar = c(margins[1], 0, 0, 0))
  par(mar = c(margins[1], 0, 0, margin_dendrogram))
  
  if (dendrogram %in% c("both", "row")) {
    # flag <- try(plot(ddr, horiz = TRUE, axes = TRUE, yaxs = "i"))
    
    flag <- try(
      ddr %>%
        dendextend::set("labels_col", colRow) %>%
        plot(yaxs = "i",horiz = TRUE))
    
    if ("try-error" %in% class(flag)) {
      cond <- attr(flag, "condition")
      if (!is.null(cond) &&
          conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?")
        stop(
          "Row dendrogram too deeply nested, recursion limit exceeded. Try increasing option(\"expressions\"=...)."
        )
    }
  } else {
    plot.new()
  }
  
  # Draw column dendrogram
  par(mar = c(margin_dendrogram, 0, if (!is.null(main)) 5 else 0, margins[2]))
  
  if (dendrogram %in% c("both", "column")) {
    # flag <- try(plot(ddc, xaxs = "i"))
    
    flag <- try(
      ddc %>%
        dendextend::set("labels_col", colCol) %>%
        plot(xaxs = "i"))
    if ("try-error" %in% class(flag)) {
      cond <- attr(flag, "condition")
      if (!is.null(cond) &&
          conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?")
        stop(
          "Column dendrogram too deeply nested, recursion limit exceeded. Try increasing option(\"expressions\"=...)."
        )
    }
  } else {
    plot.new()
  }
  
  # Add main title
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  
  # Draw color key
  if (key) {
    mar <- c(5, 4, 2, 1)
    if (!is.null(key.xlab) && is.na(key.xlab))
      mar[1] <- 2
    if (!is.null(key.ylab) && is.na(key.ylab))
      mar[2] <- 2
    if (!is.null(key.title) && is.na(key.title))
      mar[3] <- 1
    
    par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
    
    if (length(key.par) > 0)
      do.call(par, key.par)
    
    tmpbreaks <- breaks
    
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    } else {
      min.raw <- min.breaks
      max.raw <- max.breaks
    }
    
    z <- seq(min.raw, max.raw, by = min(diff(breaks) / 100))
    
    image(
      z = matrix(z, ncol = 1),
      col = col,
      breaks = tmpbreaks,
      xaxt = "n",
      yaxt = "n"
    )
    
    par(usr = c(0, 1, 0, 1))
    
    if (is.null(key.xtickfun)) {
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      xargs <- list(at = xv, labels = lv)
    } else {
      xargs <- key.xtickfun()
    }
    
    xargs$side <- 1
    do.call(axis, xargs)
    
    if (is.null(key.xlab)) {
      if (scale == "row")
        key.xlab <- "Row Z-Score"
      else if (scale == "column")
        key.xlab <- "Column Z-Score"
      else
        key.xlab <- "Value"
    }
    
    if (!is.na(key.xlab)) {
      mtext(
        side = 1,
        key.xlab,
        line = par("mgp")[1],
        padj = 0.5,
        cex = par("cex") * par("cex.lab")
      )
    }
    
    # Draw density plot if requested
    if (density.info == "density") {
      dens <- density(
        x,
        adjust = densadj,
        na.rm = TRUE,
        from = min.scale,
        to = max.scale
      )
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[!omit]
      dens$y <- dens$y[!omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      
      lines(dens$x, dens$y / max(dens$y) * 0.95, col = denscol, lwd = 1)
      
      if (is.null(key.ytickfun)) {
        yargs <- list(at = pretty(dens$y) / max(dens$y) * 0.95,
                      labels = pretty(dens$y))
      } else {
        yargs <- key.ytickfun()
      }
      
      yargs$side <- 2
      do.call(axis, yargs)
      
      if (is.null(key.title))
        key.title <- "Color Key\nand Density Plot"
      
      if (!is.na(key.title))
        title(key.title)
      
      par(cex = 0.5)
      
      if (is.null(key.ylab))
        key.ylab <- "Density"
      
      if (!is.na(key.ylab))
        mtext(
          side = 2,
          key.ylab,
          line = par("mgp")[1],
          padj = 0.5,
          cex = par("cex") * par("cex.lab")
        )
    } 
    # Draw histogram if requested
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      
      lines(hx, hy / max(hy) * 0.95, lwd = 1, type = "s", col = denscol)
      
      if (is.null(key.ytickfun)) {
        yargs <- list(at = pretty(hy) / max(hy) * 0.95, labels = pretty(hy))
      } else {
        yargs <- key.ytickfun()
      }
      
      yargs$side <- 2
      do.call(axis, yargs)
      
      if (is.null(key.title))
        key.title <- "Color Key\nand Histogram"
      
      if (!is.na(key.title))
        title(key.title)
      
      par(cex = 0.5)
      
      if (is.null(key.ylab))
        key.ylab <- "Count"
      
      if (!is.na(key.ylab))
        mtext(
          side = 2,
          key.ylab,
          line = par("mgp")[1],
          padj = 0.5,
          cex = par("cex") * par("cex.lab")
        )
    } 
    # No density info
    else {
      if (is.null(key.title))
        key.title <- "Color Key"
      
      if (!is.na(key.title))
        title(key.title)
    }
    
    # Add trace lines to key
    if (trace %in% c("both", "column")) {
      vline.vals <- scale01(vline, min.raw, max.raw)
      if (!is.null(vline)) {
        abline(v = vline.vals, col = linecol, lty = 2)
      }
    }
    
    if (trace %in% c("both", "row")) {
      hline.vals <- scale01(hline, min.raw, max.raw)
      if (!is.null(hline)) {
        abline(v = hline.vals, col = linecol, lty = 2)
      }
    }
  } else {
    # No key
    par(mar = c(0, 0, 0, 0))
    plot.new()
  }
  
  # Prepare return value
  retval$colorTable <- data.frame(
    low = retval$breaks[-length(retval$breaks)],
    high = retval$breaks[-1],
    color = retval$col
  )
  
  retval$layout <- list(
    lmat = lmat,
    lhei = lhei,
    lwid = lwid
  )
  
  # Call extra function if provided
  if (!is.null(extrafun))
    extrafun()
  
  invisible(retval)
}