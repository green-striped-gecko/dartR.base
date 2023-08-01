#' @name gl.report.locmetric
#' @title Reports summary of the slot $other$loc.metrics
#' @family matched report
#' 
#' @description
#' This function reports summary statistics (mean, minimum, average, quantiles), histograms
#' and boxplots for any loc.metric with numeric values (stored in 
#' $other$loc.metrics) to assist the decision of choosing thresholds for the filter
#' function \code{\link{gl.filter.locmetric}}.

#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param metric Name of the metric to be used for filtering [required].
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme User specified theme [default theme_dartR()].
#' @param plot.colors Vector with two color names for the borders and fill
#' [default c("#2171B5", "#6BAED6")].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]
#' @param plot.file Filename (minus extension) for the RDS plot file [Required for plot save]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].

#' @details
#'The function \code{\link{gl.filter.locmetric}} will filter out the
#'  loci with a locmetric value below a specified threshold.

#'The fields that are included in dartR, and a short description, are found
#'below. Optionally, the user can also set his/her own field by adding a vector
#' into $other$loc.metrics as shown in the example. You can check the names of
#' all available loc.metrics via: names(gl$other$loc.metrics).

#'\itemize{
#'\item SnpPosition - position (zero is position 1) in the sequence tag of the
#'defined SNP variant base.
#'\item CallRate - proportion of samples for which the genotype call is
#'non-missing (that is, not '-' ).
#'\item OneRatioRef - proportion of samples for which the genotype score is 0.
#'\item OneRatioSnp - proportion of samples for which the genotype score is 2.
#'\item FreqHomRef - proportion of samples homozygous for the Reference allele.
#'\item FreqHomSnp - proportion of samples homozygous for the Alternate (SNP)
#'allele.
#'\item FreqHets - proportion of samples which score as heterozygous, that is,
#'scored as 1.
#'\item PICRef - polymorphism information content (PIC) for the Reference allele.
#'\item PICSnp - polymorphism information content (PIC) for the SNP.
#'\item AvgPIC - average of the polymorphism information content (PIC) of the
#' reference and SNP alleles.
#'\item AvgCountRef - sum of the tag read counts for all samples, divided by the
#' number of samples with non-zero tag read counts, for the Reference allele row.
#'\item AvgCountSnp - sum of the tag read counts for all samples, divided by the
#'number of samples with non-zero tag read counts, for the Alternate (SNP) allele
#' row.
#'\item RepAvg - proportion of technical replicate assay pairs for which the
#'marker score is consistent.
#'\item rdepth - read depth.
#'}

#'\strong{ Function's output }

#' The minimum, maximum, mean and a tabulation of quantiles of the locmetric
#' values against thresholds rate are provided. Output also includes a boxplot
#' and a histogram.

#' Quantiles are partitions of a finite set of values into q subsets of (nearly)
#' equal sizes. In this function q = 20. Quantiles are useful measures because
#' they are less susceptible to long-tailed distributions and outliers.

#'  Plot colours can be set with gl.select.colors().

#'   If plot.file is specified, plots are saved to the directory specified by the user, or the global
#'   default working directory set by gl.set.wd() or to the tempdir().

#'  Examples of other themes that can be used can be consulted in:
#'   \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }

#' @author Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})

#' @examples
#' # SNP data
#' out <- gl.report.locmetric(testset.gl,metric='SnpPosition')
#' # Tag P/A data
#' out <- gl.report.locmetric(testset.gs,metric='AvgReadDepth')

#' @seealso \code{\link{gl.filter.locmetric}}

#' @export
#' @return An unaltered genlight object.

gl.report.locmetric <- function(x,
                                metric,
                                plot.display=TRUE,
                                plot.theme = theme_dartR(),
                                plot.colors = NULL,
                                plot.dir=NULL,
                                plot.file=NULL,
                                verbose = NULL) {
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
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    # check whether the field exists in the genlight object
    if (!(metric %in% colnames(x$other$loc.metrics))) {
        stop(error("  Fatal Error: name of the metric not found\n"))
    }
    if (!is.numeric(unlist(x$other$loc.metrics[metric]))) {
        stop(error("  Fatal Error: metric is not numeric\n"))
    }
    
    # DO THE JOB
    
    field <- which(colnames(x@other$loc.metrics) == metric)
    
    # get title for plots
    if (all(x@ploidy == 2)) {
        title1 <- paste0("SNP data - ", metric, " by Locus")
    } else {
        title1 <- paste0("Fragment P/A data - ", metric, " by Locus")
    }
    
    metric_df <- data.frame(x$other$loc.metrics[field])
    colnames(metric_df) <- "field"
    
    p1 <-
        ggplot(metric_df, aes(y = field)) + geom_boxplot(color = plot.colors[1], fill = plot.colors[2]) + coord_flip() + plot.theme + xlim(range = c(-1,
                                                                                                                                                     1)) + ylab(metric) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + ggtitle(title1)
    
    p2 <-
        ggplot(metric_df, aes(x = field)) + geom_histogram(bins = 50,
                                                           color = plot.colors[1],
                                                           fill = plot.colors[2]) + xlab(metric) + ylab("Count") +
        plot.theme
    
    # Print out some statistics
    stats <- summary(metric_df)
    cat("  Reporting", metric, "by Locus\n")
    cat("  No. of loci =", nLoc(x), "\n")
    cat("  No. of individuals =", nInd(x), "\n")
    cat("    Minimum      : ", stats[1], "\n")
    cat("    1st quantile : ", stats[2], "\n")
    cat("    Median       : ", stats[3], "\n")
    cat("    Mean         : ", stats[4], "\n")
    cat("    3r quantile  : ", stats[5], "\n")
    cat("    Maximum      : ", stats[6], "\n\n")
    
    # Determine the loss of loci for a given threshold using quantiles
    quantile_res <-
        quantile(metric_df$field, probs = seq(0, 1, 1 / 20),type=1)
    retained <- unlist(lapply(quantile_res, function(y) {
        res <- length(metric_df$field[metric_df$field >= y])
    }))
    pc.retained <- round(retained * 100 / nLoc(x), 1)
    filtered <- nLoc(x) - retained
    pc.filtered <- 100 - pc.retained
    df <-
        data.frame(as.numeric(sub("%", "", names(quantile_res))),
                   quantile_res,
                   retained,
                   pc.retained,
                   filtered,
                   pc.filtered)
    colnames(df) <-
        c("Quantile",
          "Threshold",
          "Retained",
          "Percent",
          "Filtered",
          "Percent")
    df <- df[order(-df$Quantile), ]
    df$Quantile <- paste0(df$Quantile, "%")
    rownames(df) <- NULL
    
    # PRINTING OUTPUTS

    # using package patchwork
      p3 <- (p1 / p2) + plot_layout(heights = c(1, 4))
      if (plot.display) {print(p3)}
      print(df)
    
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
    invisible(x)
    
}
