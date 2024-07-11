#' @name gl.filter.maf
#' @title Filters loci on the basis of minor allele frequency (MAF) or minor
#' allele count (MAC)
#' @family matched filter
#' @description
#' This script calculates the minor allele frequency for each locus and updates
#' the locus metadata for FreqHomRef, FreqHomSnp, FreqHets and MAF (if it
#' exists). It then uses the updated metadata for MAF to filter loci.
#' @param x Name of the genlight object containing the SNP data [required].
#' @param threshold Threshold MAF -- loci with a MAF less than the threshold
#' will be removed. If a value > 1 is provided it will be interpreted as MAC
#' (i.e. the minimum number of times an allele needs to be observed)
#' [default 0.01].
#' @param by.pop Whether MAF should be calculated by population [default FALSE].
#' @param pop.limit Minimum number of populations in which MAF should be less
#' than the threshold for a locus to be filtered out. Only used if by.pop
#' = TRUE. The default value is half of the populations
#' [default ceiling(nPop(x)/2)].
#' @param ind.limit Minimum number of individuals that a population should
#' contain to calculate MAF. Only used if by.pop=TRUE [default 10].
#' @param recalc Recalculate the locus metadata statistics [default FALSE].
#' @param plot.display If TRUE, histograms of base composition are displayed in
#'  the plot window [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors List of two color names for the borders and fill of the
#'  plots [default c("#2171B5", "#6BAED6")].
#' @param plot.dir Directory in which to save files
#' [default = working directory].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL].
#' @param bins Number of bins to display in histograms [default 25].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
#' @details
#' Careful consideration needs to be given to the settings to be used for this
#' function. When the filter is applied globally (i.e. \code{by.pop=FALSE}) but
#' the data include multiple population, there is the risk to remove markers
#' because the allele frequencies is low (at global level) but the allele
#' frequencies for the same markers may be high within some of the populations
#' (especially if the per-population sample size is small). Similarly, not
#' always it is a sensible choice to run this function using \code{by.pop=TRUE}
#' because allele that are rare in a population may be very common in other,
#' but the (possible) allele frequencies will depend on the sample size within
#' each population. Where the purpose of filtering for MAF is to remove possible
#'  spurious alleles (i.e. sequencing errors), it is perhaps better to filter
#'   based on the number of times an allele is observed (MAC, Minimum Allele
#'   Count), under the assumption that if an allele is observed > MAC, it is
#'   fairly rare to be an error.
#'
#' \bold{From v2.1} The threshold can take values > 1. In this case, these are
#' interpreted as a threshold for MAC.
#' @author Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' result <- gl.filter.maf(platypus.gl, threshold = 0.05, verbose = 3)
#' #result <- gl.filter.maf(platypus.gl, by.pop = TRUE, threshold = 0.05, verbose = 3)
#' @export
#' @return The reduced genlight dataset

gl.filter.maf <- function(x,
                          threshold = 0.01,
                          by.pop = FALSE,
                          pop.limit = ceiling(nPop(x) / 2),
                          ind.limit = 10,
                          recalc = FALSE,
                          plot.display = TRUE,
                          plot.theme = theme_dartR(),
                          plot.colors = NULL,
                          plot.file = NULL,
                          plot.dir = NULL,
                          bins = 25,
                          verbose = NULL) {
  hold <- x
  plot.colors.pop <- plot.colors
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # SET COLOURS
  if (is.null(plot.colors)) {
    plot.colors <- gl.select.colors(library = "brewer",
                                    palette = "Blues",
                                    select = c(7, 5), verbose = 0)
  }
  
  if (verbose == 0) {
    plot.display <- FALSE
  }
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Josh",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # Work around a bug in adegenet if genlight object is created by subsetting
  if (nLoc(x) != nrow(x@other$loc.metrics)) {
    stop(
      error(
        "The number of rows in the loc.metrics table does not match the
                number of loci in your genlight object!"
      )
    )
  }
  
  # Check for monomorphic loci
  tmp <- gl.filter.monomorphs(x, verbose = 0)
  if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {
    cat(warn("  Warning: genlight object contains monomorphic loci\n"))
  }
  
  # FUNCTION SPECIFIC ERROR CHECKING
  if (threshold >= 1)
    threshold <- threshold / (length(indNames(x)) *
                                mean(ploidy(x)))
  if (threshold > 0.5 | threshold <= 0) {
    cat(
      warn(
        "  Warning: threshold must be in the range (0,0.5], but usually
                small, set to 0.01\n"
      )
    )
    threshold <- 0.01
  }
  
  # DO THE JOB
  
  if (by.pop) {
    if (verbose >= 2) {
      cat(
        report(
          "  Removing loci with MAF <",
          threshold,
          "in at least",
          pop.limit,
          "populations and recalculating FreqHoms and FreqHets\n"
        )
      )
    }
    #x <- utils.recalc.maf(x, verbose = 0)
    pop.list <- seppop(x)
    #col=gl.select.colors(library="brewer",palette="Blues",select=c(7,5), verbose=0)
    
    # getting populations with more than ind.limit
    ind_per_pop <-
      which(unlist(lapply(pop.list, nInd)) >= ind.limit)
    pop.list <- pop.list[ind_per_pop]
    # recalculating MAF by population
    pop.list <- lapply(pop.list, utils.recalc.maf, verbose = 0)
    # getting loci with MAF < threshold
    loci.list <- lapply(pop.list, function(y) {
      y$other$loc.metrics$maf <= threshold
    })
    # getting the loci in which MAF < threshold and in at least pop.limit
    # populations
    loci.list <- Reduce("+", loci.list)
    loci.list <- which(loci.list >= pop.limit)
    
    x2 <- x[,-loci.list]
    x2@other$loc.metrics <- x@other$loc.metrics[-loci.list, ]
    
    x2 <- utils.recalc.maf(x2, verbose = 0)
  } else{
    # Recalculate the relevant loc.metrics
    if (verbose >= 2) {
      cat(
        report(
          "  Removing loci with MAF <",
          threshold,
          "over all the dataset
                and recalculating FreqHoms and FreqHets\n"
        )
      )
    }
    
    x <- utils.recalc.maf(x, verbose = 0)
    
    # Remove loci with NA count <= 1-threshold
    index <- which(x@other$loc.metrics$maf >= threshold)
    
    x2 <- x[, index]
    x2@other$loc.metrics <- x@other$loc.metrics[index, ]
    
    x2 <- utils.recalc.maf(x2, verbose = 0)
  }
  
  if (plot.display & by.pop == FALSE) {
    popn.hold <- FALSE
    maf <- NULL
    # Plot a histogram of MAF
    maf_pre <- data.frame(x@other$loc.metrics$maf)
    colnames(maf_pre) <- "maf"
    min <- min(maf_pre, threshold, na.rm = TRUE)
    min <- trunc(min * 100) / 100
    
    p1 <-
      ggplot(as.data.frame(maf_pre), aes(x = maf)) +
      geom_histogram(bins = bins,
                     color = plot.colors[1],
                     fill = plot.colors[2]) +
      coord_cartesian(xlim = c(min, 0.5)) +
      geom_vline(xintercept = threshold,
                 color = "red",
                 size = 1) +
      xlab("Pre-filter SNP MAF\nOver all populations") +
      ylab("Count") +
      plot.theme
    
    maf_post <- data.frame(x2@other$loc.metrics$maf)
    colnames(maf_post) <- "maf"
    min <- min(maf_post, threshold, na.rm = TRUE)
    min <- trunc(min * 100) / 100
    
    p2 <-
      ggplot(as.data.frame(maf_post), aes(x = maf)) +
      geom_histogram(bins = bins,
                     color = plot.colors[1],
                     fill = plot.colors[2]) +
      coord_cartesian(xlim = c(min, 0.5)) +
      geom_vline(xintercept = threshold,
                 color = "red",
                 size = 1) +
      xlab("Post-filter SNP MAF\nOver all populations") +
      ylab("Count") +
      plot.theme
  }
  
  if (plot.display & by.pop == TRUE) {
    # plots pre
    pops_maf_pre <- seppop(x)
    
    mafs_plots_pre <- lapply(pops_maf_pre, function(z) {
      pop_name <- popNames(z)
      z$other$loc.metrics <-
        as.data.frame(z$other$loc.metrics)
      z <- gl.filter.monomorphs(z, verbose = 0)
      z <- gl.recalc.metrics(z, verbose = 0)
      mafs_per_pop <- z$other$loc.metrics$maf
      
      p_temp <-
        ggplot(as.data.frame(mafs_per_pop), aes(x = mafs_per_pop)) +
        geom_histogram(bins = bins,
                       color = plot.colors[1],
                       fill = plot.colors[2]) +
        geom_vline(xintercept = threshold,
                   color = "red",
                   size = 1) +
        xlab("Pre-filter SNP MAF") +
        ylab("Count") +
        xlim(0, 0.5) +
        plot.theme +
        ggtitle(paste(popNames(z), "n =", nInd(z)))
      
      return(p_temp)
    })
    
    # plots post
    pops_maf_post <- seppop(x2)
    
    mafs_plots_post <- lapply(pops_maf_post, function(z) {
      pop_name <- popNames(z)
      z$other$loc.metrics <-
        as.data.frame(z$other$loc.metrics)
      z <- gl.filter.monomorphs(z, verbose = 0)
      z <- gl.recalc.metrics(z, verbose = 0)
      mafs_per_pop <- z$other$loc.metrics$maf
      
      p_temp <-
        ggplot(as.data.frame(mafs_per_pop), aes(x = mafs_per_pop)) +
        geom_histogram(bins = bins,
                       color = plot.colors[1],
                       fill = plot.colors[2]) +
        geom_vline(xintercept = threshold,
                   color = "red",
                   size = 1) +
        xlab("Post-filter SNP MAF\n") +
        ylab("Count") +
        xlim(0, 0.5) +
        plot.theme
      
      return(p_temp)
    })
    
    plots_pops_merge <-
      lapply(1:length(pops_maf_pre), function(y) {
        plot_temp <- mafs_plots_pre[[y]] / mafs_plots_post[[y]]
        return(plot_temp)
      })
    
    # Check for status -- any populations with ind > ind.limit; and is
    #nPop > 1
    ind_per_pop <- unlist(lapply(pops_maf_pre, nInd))
    
    test_pop <-
      as.data.frame(cbind(pop = names(ind_per_pop), ind_per_pop))
    test_pop$ind_per_pop <- as.numeric(test_pop$ind_per_pop)
    
    maf_pre <- data.frame(x@other$loc.metrics$maf)
    colnames(maf_pre) <- "maf"
    
    maf_post <- data.frame(x2@other$loc.metrics$maf)
    colnames(maf_post) <- "maf"
    
    # testing which populations comply with thresholds
    popn.hold <-
      test_pop[which(test_pop$ind_per_pop >= ind.limit), "pop"]
    
    names(plots_pops_merge) <- popn.hold
    
    mafs_plots_print <- plots_pops_merge[popn.hold]
    
    if (length(popn.hold) > 1) {
      p_all_pre <-
        ggplot(as.data.frame(maf_pre), aes(x = maf)) +
        geom_histogram(bins = bins,
                       color = plot.colors[1],
                       fill = plot.colors[2]) +
        geom_vline(xintercept = threshold,
                   color = "red",
                   size = 1) +
        xlab("Pre-filter SNP MAF\nOver all populations") +
        ylab("Count") +
        xlim(0, 0.5) +
        plot.theme
      
      p_all_post <-
        ggplot(as.data.frame(maf_post), aes(x = maf)) +
        geom_histogram(bins = bins,
                       color = plot.colors[1],
                       fill = plot.colors[2]) +
        geom_vline(xintercept = threshold,
                   color = "red",
                   size = 1) +
        xlab("Post-filter SNP MAF\nOver all populations") +
        ylab("Count") +
        xlim(0, 0.5) +
        plot.theme
      
      p2 <- p_all_pre / p_all_post
      p3 <- mafs_plots_print
    }
    
    if (length(popn.hold) == 0) {
      if (verbose >= 1) {
        cat(
          important(
            "  No populations met minimum limits on number of
                        individuals or loci, reporting for overall\n"
          )
        )
      }
      p_all_pre <-
        ggplot(as.data.frame(maf_pre), aes(x = maf)) +
        geom_histogram(bins = bins,
                       color = plot.colors[1],
                       fill = plot.colors[2]) +
        geom_vline(xintercept = threshold,
                   color = "red",
                   size = 1) +
        xlab("Pre-filter SNP MAF\nOver all populations") +
        ylab("Count") +
        xlim(0, 0.5) +
        plot.theme
      
      p_all_post <-
        ggplot(as.data.frame(maf_post), aes(x = maf)) +
        geom_histogram(bins = bins,
                       color = plot.colors[1],
                       fill = plot.colors[2]) +
        geom_vline(xintercept = threshold,
                   color = "red",
                   size = 1) +
        xlab("Post-filter SNP MAF\nOver all populations") +
        ylab("Count") +
        xlim(0, 0.5) +
        plot.theme
      
      p3 <- p_all_pre / p_all_post
    }
    
    if (length(popn.hold) == 1) {
      if (verbose >= 3) {
        cat(
          important(
            "  Only one population met minimum limits on number of
                        individuals or loci\n"
          )
        )
      }
      
      p_all_pre <-
        ggplot(as.data.frame(maf_pre), aes(x = maf)) +
        geom_histogram(bins = bins,
                       color = plot.colors[1],
                       fill = plot.colors[2]) +
        geom_vline(xintercept = threshold,
                   color = "red",
                   size = 1) +
        xlab("Pre-filter SNP MAF") +
        ylab("Count") +
        xlim(0, 0.5) +
        plot.theme
      
      p_all_post <-
        ggplot(as.data.frame(maf_post), aes(x = maf)) +
        geom_histogram(bins = bins,
                       color = plot.colors[1],
                       fill = plot.colors[2]) +
        geom_vline(xintercept = threshold,
                   color = "red",
                   size = 1) +
        xlab("Post-filter SNP MAF") +
        ylab("Count") +
        xlim(0, 0.5) +
        plot.theme
      
      p3 <- p_all_pre / p_all_post
      
    }
    
    if (nPop(x2) == 1) {
      if (verbose >= 1) {
        cat(important("  Only one population specified\n"))
      }
      title.str <-
        paste("Minor Allele Frequency\n", pop(x2)[1])
      
      p_all_pre <-
        ggplot(as.data.frame(maf_pre), aes(x = maf)) +
        geom_histogram(bins = bins,
                       color = plot.colors[1],
                       fill = plot.colors[2]) +
        geom_vline(xintercept = threshold,
                   color = "red",
                   size = 1) +
        xlab("Pre-filter SNP MAF") +
        ylab("Count") +
        xlim(0, 0.5) +
        plot.theme
      
      p_all_post <-
        ggplot(as.data.frame(maf_post), aes(x = maf)) +
        geom_histogram(bins = bins,
                       color = plot.colors[1],
                       fill = plot.colors[2]) +
        geom_vline(xintercept = threshold,
                   color = "red",
                   size = 1) +
        xlab("Post-filter SNP MAF") +
        ylab("Count") +
        xlim(0, 0.5) +
        plot.theme
      
      p3 <- p_all_pre / p_all_post
    }
  }
  
  if (recalc) {
    # Recalculate all metrics(flags reset in utils scripts)
    x2 <- gl.recalc.metrics(x2, verbose = verbose)
  } else {
    # Reset the flags as FALSE for all metrics except MAF (dealt with elsewhere)
    x2@other$loc.metrics.flags$AvgPIC <- FALSE
    x2@other$loc.metrics.flags$OneRatioRef <- FALSE
    x2@other$loc.metrics.flags$OneRatioSnp <- FALSE
    x2@other$loc.metrics.flags$PICRef <- FALSE
    x2@other$loc.metrics.flags$PICSnp <- FALSE
    x2@other$loc.metrics.flags$FreqHets <- FALSE
    x2@other$loc.metrics.flags$FreqHomRef <- FALSE
    x2@other$loc.metrics.flags$FreqHomSnp <- FALSE
    x2@other$loc.metrics.flags$CallRate <- FALSE
    x2@other$loc.metrics.flags$allna <- FALSE
  }
  
  # REPORT A SUMMARY
  if (verbose >= 3) {
    cat("  Summary of filtered dataset\n")
    cat("  MAF for loci >", threshold, "\n")
    cat("  Initial number of loci:", nLoc(x), "\n")
    cat("  Number of loci deleted:", nLoc(x) - nLoc(x2), "\n")
    cat("  Final number of loci:", nLoc(x2), "\n")
  }
  
  # PRINTING OUTPUTS using package patchwork
  if (plot.display) {
    if (length(popn.hold) > 1 & by.pop == TRUE) {
      suppressWarnings(print(p2))
      suppressWarnings(print(p3))
      p3 <- c(p3, p2)
    } else{
      p3 <- (p1 / p2) + plot_layout(heights = c(1, 1))
      suppressWarnings(print(p3))
    }
  }
  
  # Optionally save the plot ---------------------
  
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p3,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }
  
  # ADD TO HISTORY
  nh <- length(x2@other$history)
  x2@other$history[[nh + 1]] <- match.call()
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(x2)
}
