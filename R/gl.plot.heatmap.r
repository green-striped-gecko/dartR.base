#' @name gl.plot.heatmap
#' @title Represents a distance matrix as a heatmap
#' @family graphics

#' @description
#' The script plots a heat map to represent the distances in the distance or
#' dissimilarity matrix. This function is a wrapper for
#' \link[gplots]{heatmap.2} (package gplots).
#' @param D Name of the distance matrix or class fd object [required].
#' @param x Genlight object to extract population information [default NULL].
#' @param palette.divergent A divergent palette for the distance values
#'  [default gl.colors("div")].
#' @param palette_discrete The color of populations [default NULL].
#' @param dendrogram Character string indicating whether to draw 'none',
#' 'row', 'column' or 'both' dendrograms [default "column"].
#' @param plot.out A boolean that indicates whether to plot the results [default TRUE].
#' @param legend.print Whether to create legend (only if x is provided) 
#' [default TRUE].
#' @param legendx x coordinates for the legend[default 0].
#' @param legendy y coordinates for the legend[default 1].
#' @param label.size Specify the size of the population labels [default 0.75].
#' @param legend.title Legend title [default "Populations"].
#' @param diag.na Logical. If TRUE, the diagonal elements of the distance matrix are set to NA [default FALSE].
#' @param margins Numeric vector of length 2 containing the margins for column
#' and row names, respectively [Default = c(10, 10)].
#' @param na.color Color to use for missing value (NA) [default "grey"].
#' @param revC Reverse column order. Logical value.
#' @param symbreaks Symmetric breaks setting. Logical value. Default: FALSE
#' @param trace Draw trace line. Must be one of: "column", "row", "both", "none". Default: "none"
#' @param tracecol Trace line colour. Default: "cyan"
#' @param cexRow Row axis label scale. Integer value.
#' @param cexCol Column axis label scale. Integer value.
#' @param srtRow Row labels' rotation angle. Integer value.
#' @param srtCol Column labels' rotation angle. Integer value. Default: 90
#' @param offsetRow Row labels' offset. Integer value. Default: 0.5
#' @param offsetCol Column labels' offset. Integer value. Default: 0.5
#' @param key Show colour-key. Logical value. Default: TRUE
#' @param keysize Size of colour key. Numeric value (>= 0). Default: 1.5
#' @param density.info Density plot for colour key. Must be one of: "histogram", "density", "none". Default: "none"
#' @param denscol Density plot colour. Default: "cyan"
#' @param symkey Symmetric colour key. Logical value. Default: FALSE
#' @param densadj Density plot scaling. Numeric value between 0 and 10. Default: 0.25
#' @param key.title Title for colour key. Character string.
#' @param key.xlab X title for colour key. Character string.
#' @param key.ylab Y title for colour key. Character string.
#' @param main Main plot title. Character string.
#' @param xlab X-axis label. Character string.
#' @param ylab Y-axis label. Character string.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity]
#' @param ... Parameters passed to function
#' \link[gplots]{heatmap.2} (package gplots)
#'
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#'    gl <- testset.gl[1:10,]
#'    if (isTRUE(getOption("dartR_fbm"))) gl <- gl.gen2fbm(gl)
#'    D <- dist(as.matrix(gl),upper=TRUE,diag=TRUE)
#'    gl.plot.heatmap(D)
#'    if (isTRUE(getOption("dartR_fbm"))) possums.gl <- gl.gen2fbm(possums.gl)
#'    D2 <- gl.dist.pop(possums.gl)
#'    gl.plot.heatmap(D2)
#'    if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#'    D3 <- gl.fixed.diff(testset.gl)
#'    gl.plot.heatmap(D3)
#'    }
#'    if ((requireNamespace("gplots", quietly = TRUE))) {
#'    D2 <- gl.dist.pop(possums.gl)
#'    gl.plot.heatmap(D2)
#'    }
#' @importFrom graphics legend
#' @importFrom gtools invalid 
#' @export
#' @return returns no value (i.e. NULL)

gl.plot.heatmap <- function(D,
                            x = NULL,
                            palette.divergent = gl.colors("div"),
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
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  # DO THE JOB
  
  colors_pops <- NULL
  legend_text <-  ""
  legend_color <-  "white"
  
  if(is.null(x)){
    legend.print <- FALSE
  }
  
  if (!is.null(x)) {
    # assigning colors to populations
    if (!is.null(palette_discrete)) {
      # if pop colors is a palette
      if (is(palette_discrete, "function")) {
        colors_pops <- palette_discrete(nPop(x))
      }
      # if pop colors is a vector
      if (!is(palette_discrete, "function")) {
        colors_pops <- palette_discrete
      }
    } else{
      colors_pops <- gl.select.colors(x, verbose = 0)
    }
    
    names(colors_pops) <- as.character(popNames(x))
    df_colors_temp_1 <- data.frame(ind = indNames(x), pop = as.character(pop(x)))
    
    df_colors_temp_2 <- data.frame(pop = names(colors_pops), color = colors_pops)
    df_colors <- merge(df_colors_temp_1, df_colors_temp_2, by = "pop")
    
    df_mat <- data.frame(ind = colnames(as.matrix(D)), order = 1:ncol((as.matrix(D))))
    df_colors_2 <- merge(df_mat, df_colors, by = "ind")
    
    df_colors_2 <- df_colors_2[order(df_colors_2$order), ]
    
    colors_pops <- df_colors_2$color
    
    legend_text <-  unique(df_colors_2$pop)
    legend_color <-  unique(df_colors_2$color)
    
  }
  
  if (datatype == "dist" | datatype == "matrix") {
    D <- as.dist(D)
    m <- as.matrix(D)
    
    if (!is.null(x) && ncol(m) != nInd(x)) {
      colors_pops <- NULL
      legend.print <- FALSE
    }
    
    if (diag.na) diag(m) <- NA
    
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
      par(mar = c(1, 1, 1, 1))
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
  }
  
  if (datatype == "fd") {
    m <- as.matrix(D$fd)
    
    if (!is.null(x) && ncol(m) != nInd(x)) {
      colors_pops <- NULL
    }
    
    if (diag.na) diag(m) <- NA
    
    if (plot.out) {
    p3 <-
      utils.heatmap(
        m,
        col = palette.divergent(255),
        dendrogram = dendrogram,
        margins = margins,
        na.color = na.color,
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
    }
  }
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  if (plot.out) {
  invisible(p3)
}
}
