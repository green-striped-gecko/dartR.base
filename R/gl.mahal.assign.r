#' @name gl.mahal.assign
#' @title Assigns individuals to populations with an associated probability
#' @family matched reports
#'
#' @description
#' Uses Mahalanobis Distances between individuals and group centroids to
#' calculate probability of group membership
#'
#' @param x Name of the genlight object [required].
#' @param plot.display If TRUE, resultant plots are displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors List of two color names for the borders and fill of the
#'  plots [default c("#2171B5","#6BAED6")].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#'
#' @details 
#' The function generates 200 simulated individuals for each group (=population in the
#' genlight object) drawing from the observed allele frequencies for each group.
#' The group centroids and covariance matricies are calculated for these simulated
#' groups. The covariance matrix is inverted using package MASS::ginv to overcome the
#' singularities that would otherwise arise with typical SNP data. Mahanobilis
#' Distances are calculated using stats::mahanalobis for each individual in the
#' dataset and associated Chi Square probabilities of group membership are calculated for
#' each individual in the original genlight object. The resultant table can be
#' used for decisions on group membership. A special group (=population in the
#' genlight object) called 'unknowns' can be used to specifically identify 
#' individuals with unknown group membership.
#'
#'  A color vector can be obtained with gl.select.colors() and then passed to the function
#'  with the plot.colors parameter.

#' Themes can be obtained from in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'   If a plot.file is given, the ggplot arising from this function is saved as an "RDS" 
#' binary file using saveRDS(); can be reloaded with readRDS(). A file name must be 
#' specified for the plot to be saved.

#'  If a plot directory (plot.dir) is specified, the ggplot binary is saved to that
#'  directory; otherwise to the tempdir(). 

#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
# @examples
#'
#' @importFrom MASS ginv
#' @importFrom stats cov mahalanobis
#' @export
#' @return The unchanged genlight object
#' 

gl.mahal.assign <- function(x,
                            plot.display=TRUE,
                            plot.theme = theme_dartR(),
                            plot.colors = NULL,
                            plot.file=NULL,
                            plot.dir=NULL,
                            verbose = NULL) {

# PRELIMINARIES -- checking 
  
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
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = c("genlight", "SNP", "SilicoDArT"), verbose = verbose)

# DO THE JOB
  
# Reduce the dataset to a dense matrix

x <- gl.filter.allna(x, verbose=0)
x <- gl.filter.monomorphs(x,verbose=0)
x <- gl.filter.callrate(x,threshold=0.95,verbose=0)
x <- gl.impute(x,verbose=0)
mat <- as.matrix(x)

# For each population

for (i in 1:nPop(x)){

  pop.name <- popNames(x)[i]
  gl <- gl.keep.pop(x,pop.list=pop.name,verbose=0)
  
  # Generate a random dataset

  gl.expanded <- gl.sim.genotypes(gl, n.ind=min(nLoc(gl),200),verbose=0)

# Compute means and covariances
  mat.expanded <- as.matrix(gl.expanded)
  means <- colMeans(mat.expanded)
  covmat <- stats::cov(mat.expanded)

# Invert the covariance matrix using MASS (mahalanobis algorithm too fussy)

  inv <- MASS::ginv(covmat)

# Calculate the Mahal distance for each point in the original genlight object

  tmp <- stats::mahalanobis(mat, center=means, cov=inv, inverted=TRUE)
  tmp <- as.data.frame(tmp)
  tmp$id <- row.names(tmp)
  tmp$prob <- round(pchisq(tmp$tmp, df=59, lower.tail=FALSE),4)
  tmp$tmp <- NULL
  #tmp <- tmp[,c(2,1)]

# Progressive data addition  
  
  if(i>1){
    result <- merge(result,tmp,by="id")
  } else {
    result <- tmp
  }

}
names(result) <- c("id",popNames(x))

# FLAG SCRIPT END ---------------

if (verbose >= 1) {
  cat(report("Completed:", funname, "\n"))
}

return(result)

}

