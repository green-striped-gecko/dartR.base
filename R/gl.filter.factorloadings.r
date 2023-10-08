#' @name gl.filter.factorloadings
#' @title Filters loci based on factor loadings for a PCA or PCoA
#' @family matched filters

#' @description
#' Extracts the factor loadings from a glPCA object (generated by gl.pcoa) and 
#' filters loci based on a user specified threshold for the ABSOLUTE value of the 
#' factor loadings.
#'  
#' @param x Name of the genlight object containing the SNP data or the 
#' SilocoDArT data [required].
#' @param pca Name of the glPCA object containing factor loadings [required].
#' @param axis Axis in the ordination used to display the factor loadings [default 1]
#' @param threshold Remove [retain=FALSE] or retain [retain=TRUE] only those loci that
#' load high (greater than the threshold) against the specified axis [required].
#' @param retain If true, the resultant genlight object holds only the loci that load
#' high on the specified axis; if FALSE, the resultant genlight object has the
#' loci loading high on the specified axis filtered out [default FALSE].
#' @param plot.display If TRUE, resultant plots are displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors List of two color names for the borders and fill of the
#'  plots [default c("#2171B5","#6BAED6")].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param bins Number of bins to display in histograms [default 25].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#' @param ... Parameters passed to function \link[ggplot2]{ggsave}, 
#'  such as width and height, when the ggplot is to be saved.

#' @details 
#' The function extracts the factor loadings for a given axis from a PCA object 
#' generated by gl.pcoa and then filters loci on the basis of a user specified
#' threshold. The threshold value is decided using gl.report.factorloadings.  The function
#' can be used to filter out loci that load high with a particular axis or alternatively
#' if retain=TRUE, to retain loci that load high on a specified axis.  
#' 
#'  A color vector can be obtained with gl.select.colors() and then passed to the function
#'  with the plot.colors parameter.
#'  
#' Themes can be obtained from in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'   If a plot.file is given, the ggplot arising from this function is saved as an "RDS" 
#' binary file using saveRDS(); can be reloaded with readRDS(). A file name must be 
#' specified for the plot to be saved.

#'  If a plot directory (plot.dir) is specified, the ggplot binary is saved to that
#'  directory; otherwise to the tempdir(). 
#'  
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' pca <- gl.pcoa(testset.gl)
#' gl.report.factorloadings(pca = pca)
#' gl2 <- gl.filter.factorloadings(pca=pca,x=testset.gl,threshold=0.2)
#' 
#' @export
#' @return The unchanged genlight object
#' 
gl.filter.factorloadings <- function(x,
                                     pca,
                                     axis=1,
                                     threshold,
                                     retain=FALSE,
                                     plot.display=TRUE,
                                     plot.theme = theme_dartR(),
                                     plot.colors = NULL,
                                     plot.file=NULL,
                                     plot.dir=NULL,
                                     bins=25,
                                     verbose = NULL,
                                     ...) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  if(verbose==0){plot.display <- FALSE}
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir,verbose=0)
  
  # SET COLOURS
  if(is.null(plot.colors)){
    plot.colors <- c("#2171B5", "#6BAED6")
  } else {
    if(length(plot.colors) > 2){
      if(verbose >= 2){cat(warn("  More than 2 colors specified, only the first 2 are used\n"))}
      plot.colors <- plot.colors[1:2]
    }
  }
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)
  
  # CHECK DATATYPE FOR THE genlight OBJECT
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # CHECK DATATYPE FOR THE pca OBJECT
  datatype <- class(pca)
  if(datatype != "glPca"){
    cat(error("To report factor loadings, require a glPca object\n"))
    stop()
  } else {
    if(verbose >= 2){cat(report("  Reading a glPca object\n"))}
  }
  
  # DO THE JOB
  
  # Remove monomorphic loci from the genlight object, because gl.pcoa does
  # this when undertaking the orination
  x <- gl.filter.monomorphs(x,verbose=0)
  
  # Pull the factor loadings into a dataframe
  factor.loadings <- data.frame(pca$loadings[,axis])
  rownames(factor.loadings) <- locNames(x)
  df <- cbind(rownames(factor.loadings),factor.loadings[,1])
  df <- data.frame(df)
  colnames(df) <- c("locus","loading")
  df$loading <- as.numeric(df$loading)
  if(retain){
    tmp <- df[abs(df$loading) >= threshold,]
    x2 <- gl.keep.loc(x,loclist<-tmp$locus,verbose=0)
    if(verbose >= 2){cat(report("  Retaining",nLoc(x2),"of",nLoc(x),"loci with loadings greater than or equal to",threshold,"\n"))}
  } else {
    tmp <- df[abs(df$loading) < threshold,]
    x2 <- gl.keep.loc(x,loclist<-tmp$locus,verbose=0)
    if(verbose >= 2){cat(report("  Retaining",nLoc(x2),"of",nLoc(x),"loci with loadings less than",threshold,"\n"))}
  }
  
  # Prepare the plots
  # get title for plots
  xlabel1 <- paste("Distribution of factor loadings for axis PRE-FILTER",axis)
  xlabel2 <- paste("Distribution of factor loadings for axis POST-FILTER",axis)
  
  # Calculate minimum and maximum graph cutoffs for callrate
  min <- min(df$loading)
  max <- max(df$loading)
  
  p1 <-
    ggplot(df, aes(x = loading)) +
    geom_histogram(bins = bins,
                   color = plot.colors[1],
                   fill = plot.colors[2]) +
    coord_cartesian(xlim = c(min, max)) + 
    geom_vline(xintercept = threshold,
               color = "red",
               linewidth = 1) + 
    geom_vline(xintercept = -threshold,
               color = "red",
               linewidth = 1) +
    xlab(xlabel1) + 
    ylab("Count") +
    plot.theme
  
  p2 <-
    ggplot(tmp, aes(x = loading)) + 
    geom_histogram(bins = bins,
                   color = plot.colors[1],
                   fill = plot.colors[2]) +
    coord_cartesian(xlim = c(min, max)) +
    # geom_vline(xintercept = threshold,
    #            color = "red",
    #            linewidth = 1) +
    xlab(xlabel2) +
    ylab("Count") +
    plot.theme
  
  # PRINTING OUTPUTS using package patchwork
  p3 <- (p1 / p2) + plot_layout(heights = c(1, 1))
  if (plot.display) {print(p3)}
  
  if(!is.null(plot.file)){
    tmp <- utils.plot.save(p3,
                           dir=plot.dir,
                           file=plot.file,
                           verbose=verbose)
  }
  
  # FLAG SCRIPT END 
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  invisible(x2)
}