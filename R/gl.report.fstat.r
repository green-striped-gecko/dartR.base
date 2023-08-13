#' @name gl.report.fstat
#' @title
#' Reports various statistics of genetic differentiation between
#' populations with confident intervals
#' @family matched reports
#' @description
#' This function calculates four genetic differentiation between populations
#' statistics (see the "Details" section for further information).
#'
#' \itemize{
#' \item "Fst" - Nei's Gst (Nei, 1987, pp. 164–165).
#' \item "Fstp" - Fst corrected for sample size (Nei, 1987, pp. 164–165).
#' \item "Dest" - Jost’s D (Jost, 2008).
#' \item "Gst_H" - Gst standardized by the maximum level that it can obtain for
#' the observed amount of genetic variation (Hedrick 2005).
#' }
#'
#' Confident Intervals are obtained using bootstrapping.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param nboots Number of bootstrap replicates to obtain confident intervals
#' [default 0].
#' @param conf The confidence level of the required interval  [default 0.95].
#' @param CI.type Method to estimate the confident intervals (CI). One of
#' "norm", "basic", "perc" or "bca" [default "bca"].
#' @param parallel The type of parallel operation to be used. One of "no",
#'  "multicore" or "snow". See details in the "Parallel operation" section from
#'  the function \link[boot]{boot} (package boot) [default "no"].
#' @param ncpus Number of processes to be used in parallel operation [default 1].
#' @param plot.stat Statistic to plot. One of "Fst","Fstp","Dest" or "Gst_H"
#' [default "Fstp"].
#' @param plot.display If TRUE, a heatmap of the pairwise static chosen is
#'  displayed in the plot window [default TRUE].
#' @param palette.divergent A color palette function for the heatmap plot
#'  [default gl.colors("div")].
#' @param font.size Size of font for the labels of horizontal and vertical axes
#' of the heatmap [default 0.5].
#' @param plot.dir Directory in which to save files [default = working directory].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#' @param ... Parameters passed to function \link[gplots]{heatmap.2} (package
#' gplots).
#' @details
#'
#'  Even though Fst and its relatives can predict evolutionary processes
#'  (Holsinger & Weir, 2009), they are not true measures of genetic
#'   differentiation in the sense that they are dependent on the diversity
#'    within populations (Meirmans & Hedrick, 2011), the number of populations
#'    analysed (Alcala & Rosenberg, 2017) and are not monotonic
#'    (Sherwin et al., 2017). Recent approaches have been developed to
#'    accommodate these mathematical restrictions (G'ST; "Gst_H"; Hedrick, 2005,
#' and Jost's D; "Dest"; Jost, 2008). More recently, novel approaches based on
#' information theory (Mutual Information; Sherwin et al., 2017) and allele
#' frequencies (Allele Frequency Difference; Berner, 2019) have distinct
#' properties that make them valuable resources to interpret genetic
#'  differentiation between populations.
#'
#'     Note that each measure of genetic differentiation has advantages and
#'     drawbacks, and the decision of using a particular measure is usually
#'     based on the research question.
#'
#'     \strong{Statistics calculated}
#'
#'     The equations to calculate the statistics are shown below.
#'
#'      \itemize{
#'      \item
#'     "Ho" - Observed heterozygosity corrected for sample size
#'     (Nei, 1987, pp. 164–165) is calculated as:
#'
#'     \figure{Ho-equation.jpg}
#'
#'     where Pkii represents the proportion of homozygote i in sample k and np
#'     the number of samples (i.e. populations).
#'
#'     \item
#'     "Hs" - Expected heterozygosity corrected for sample size
#'     (Nei, 1987, pp. 164–165) is calculated as:
#'
#'     \figure{He-equation.jpg}
#'
#'     \figure{He-equation-2.jpg}
#'
#'     \item
#'     "Ht" - Overall heterozygosity corrected for sample size
#'     (Nei, 1987, pp. 164–165) is calculated as:
#'
#'     \figure{Ht-equation.jpg}
#'
#'     \figure{Ht-equation-2.jpg}
#'
#'     \item
#'     "Dst" - Amount of heterozygosity among samples (Nei, 1987, pp. 164–165)
#'     is calculated as:
#'
#'     \figure{Dst-equation.jpg}
#'
#'     \item
#'     "Htp" - Overall heterozygosity corrected for sample size
#'     (Nei, 1987, pp. 164–165) is calculated as:
#'
#'     \figure{Htp-equation.jpg}
#'
#'     \item
#'     "Dstp" - Amount of heterozygosity among samples corrected for sample size
#'      (Nei, 1987, pp. 164–165) is calculated as:
#'
#'     \figure{Dstp-equation.jpg}
#'
#'     \item
#'     "Fst" - Nei's Gst (Nei, 1987, pp. 164–165) is calculated as:
#'
#'     \figure{Fst-equation.jpg}
#'
#'     \item
#'     "Fstp" - Fst corrected for sample size (Nei, 1987, pp. 164–165) is
#'     calculated as:
#'
#'     \figure{Fstp-equation.jpg}
#'
#'     \item
#'     "Fis" - Inbreeding coefficient is calculated as:
#'
#'     \figure{Fis-equation.jpg}
#'
#'     \item
#'     "Dest" - Jost’s D (Jost, 2008) is calculated as:
#'
#'     \figure{Dest-equation.jpg}
#'
#'     \item
#'     "Gst_max" - The maximum level that Gst can obtain for the observed amount
#'      of genetic variation (Hedrick 2005) is calculated as:
#'
#'     \figure{GstMax-equation.jpg}
#'
#'     where k is the number of subpopulations.
#'
#'     \item
#'     "Gst_H" - Gst standardized by the maximum level that it can obtain for the
#'     observed amount of genetic variation (Hedrick 2005) is calculated as:
#'
#'     \figure{Gst-H.jpg}
#'
#'     }
#'
#'  \strong{Confident Intervals}
#'
#' The uncertainty of a parameter, in this case the mean of the statistic, can
#' be summarised by a confidence interval (CI) which includes the true parameter
#' value with a specified probability (i.e. confidence level; the parameter
#' "conf" in this function).
#'
#' In this function, CI are obtained using Bootstrap which is an inference
#' method that samples with replacement the data (i.e. loci) and calculates the
#'  statistics every time.
#'
#'  This function uses the function \link[boot]{boot} (package boot) to perform
#'  the bootstrap replicates and the function \link[boot]{boot.ci}
#'  (package boot) to perform the calculations for the CI.
#'
#'  Four different types of nonparametric CI can be calculated
#'   (parameter "CI.type" in this function):
#'   \itemize{
#'    \item First order normal approximation interval ("norm").
#'    \item Basic bootstrap interval ("basic").
#'    \item Bootstrap percentile interval ("perc").
#'    \item Adjusted bootstrap percentile interval ("bca").
#'    }
#'
#' The studentized bootstrap interval ("stud") was not included in the CI types
#'  because it is computationally intensive, it may produce estimates outside
#'  the range of plausible values and it has been found to be erratic in
#'  practice, see for example the "Studentized (t) Intervals" section in:
#'
#'    \url{https://www.r-bloggers.com/2019/09/understanding-bootstrap-confidence-interval-output-from-the-r-boot-package/}
#'
#'     Nice tutorials about the different types of CI can be found in:
#'
#'     \url{https://www.datacamp.com/tutorial/bootstrap-r}
#'
#'     and
#'
#'    \url{https://www.r-bloggers.com/2019/09/understanding-bootstrap-confidence-interval-output-from-the-r-boot-package/}
#'
#'      Efron and Tibshirani (1993, p. 162) and Davison and Hinkley
#'      (1997, p. 194) suggest that the number of bootstrap replicates should
#'      be between 1000 and 2000.
#'
#'  \strong{It is important} to note that unreliable confident intervals will be
#'   obtained if too few number of bootstrap replicates are used.
#'   Therefore, the function \link[boot]{boot.ci} will throw warnings and errors
#'    if bootstrap replicates are too few. Consider increasing then number of
#'    bootstrap replicates to at least 200.
#'
#'    The "bca" interval is often cited as the best for theoretical reasons,
#'    however it may produce unstable results if the bootstrap distribution
#'     is skewed or has extreme values. For example, you might get the warning
#'     "extreme order statistics used as endpoints" or the error "estimated 
#'     adjustment 'a' is NA". In this case, you may want to use more bootstrap 
#'     replicates or a different method or check your data for outliers.
#'
#'    The error "estimated adjustment 'w' is infinite" means that the estimated
#'    adjustment ‘w’ for the "bca" interval is infinite, which can happen when
#'    the empirical influence values are zero or very close to zero. This can
#'    be caused by various reasons, such as:
#'
#'    The number of bootstrap replicates is too small, the statistic of interest
#'     is constant or nearly constant across the bootstrap samples, the data
#'     contains outliers or extreme values.
#'
#'     You can try some possible solutions, such as:
#'
#' Increasing the number of bootstrap replicates, using a different type of
#' bootstrap confidence interval or removing or transforming the outliers or
#'  extreme values.
#'
#'  \strong{Plotting}
#'
#'  The plot can be customised by including any parameter(s) from the function
#'  \link[gplots]{heatmap.2} (package gplots).
#'
#'  For the color palette you could try for example:
#'
#'  > \code{library(viridis)}
#'
#'  > \code{res <- gl.report.fstat(platypus.gl, palette.divergent = viridis)}
#'
#'   If a plot.file is given, the ggplot arising from this function is saved as an "RDS" #' binary file using saveRDS(); can be reloaded with readRDS(). A file name must be 
#' specified for the plot to be saved.
#'  If a plot directory (plot.dir) is specified, the ggplot binary is saved to that
#'  directory; otherwise to the tempdir(). 
#'
#'  Your plot might not shown in full because your 'Plots' pane is too small
#'  (in RStudio).
#'  Increase the size of the 'Plots' pane before running the function.
#'  Alternatively, use the parameter 'plot.file' to save the plot to a file.
#'
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' res <- gl.report.fstat(platypus.gl)
#'
#' @references
#' \itemize{
#' \item
#' Alcala, N., & Rosenberg, N. A. (2017). Mathematical constraints on FST:
#' Biallelic markers in arbitrarily many populations. Genetics (206), 1581-1600.
#' \item
#' Berner, D. (2019). Allele frequency difference AFD–an intuitive alternative
#' to FST for quantifying genetic population differentiation. Genes, 10(4), 308.
#' \item
#' Davison AC, Hinkley DV (1997). Bootstrap Methods and their Application.
#'  Cambridge University Press: Cambridge.
#' \item
#' Efron, B. (1979). Bootstrap methods: Another look at the jackknife. Annals of
#' Statistics 7, 1–26.
#' \item
#' Efron B, Tibshirani RJ (1993). An Introduction to the Bootstrap. Chapman and
#'  Hall: London.
#' \item
#' Hedrick, P. W. (2005). A standardized genetic differentiation measure.
#' Evolution, 59(8), 1633-1638.
#' \item
#' Holsinger, K. E., & Weir, B. S. (2009). Genetics in geographically structured
#'  populations: defining, estimating and interpreting FST. Nature Reviews
#'  Genetics, 10(9), 639- 650.
#'  \item
#'  Jost, L. (2008). GST and its relatives do not measure differentiation.
#'  Molecular Ecology, 17(18), 4015-4026.
#'  \item
#'  Meirmans, P. G., & Hedrick, P. W. (2011). Assessing population structure:
#'  FST and related measures. Molecular Ecology Resources, 11(1), 5-18.
#'  \item
#'  Nei, M. (1987). Molecular evolutionary genetics: Columbia University Press.
#'  \item
#'  Sherwin, W. B., Chao, A., Jost, L., & Smouse, P. E. (2017). Information
#'  theory broadens the spectrum of molecular ecology and evolution. Trends in
#'   Ecology & Evolution, 32(12), 948-963.
#' }
#' @export
#' @return Two lists, the first list contains matrices with genetic statistics
#' taken pairwise by population, the second list contains tables with the
#' genetic statistics for each pair of populations. If nboots > 0, tables with
#' the four statistics calculated with Low Confidence Intervals (LCI) and High
#' Confidence Intervals (HCI).
#'
# ----------------------
# Function

gl.report.fstat <- function(x,
                            nboots = 0,
                            conf = 0.95,
                            CI.type = "bca",
                            parallel = "no",
                            ncpus = 1,
                            plot.stat = "Fstp",
                            plot.display = TRUE,
                            palette.divergent = gl.colors("div"),
                            font.size = 0.5,
                            plot.dir = NULL,
                            plot.file = NULL,
                            verbose = NULL,
                            ...) {
  # PRELIMINARIES -- checking ----------------
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # keeping populations with more than 1 individuals
  pop_names <- popNames(x)[which(table(pop(x)) > 1)]
  
  if (length(pop_names) < nPop(x)) {
    if (verbose >= 2) {
      cat(warn("   Keeping populations with more than one individuals.\n"))
    }
    x <- gl.keep.pop(x,
                     pop.list = pop_names,
                     verbose = verbose)
  }
  
  #converting to dartR object
  class(x) <- "dartR"
  
  # bootstrapping function
  pop.diff <- function(x, indices) {
    pop.diff_fun <- function(y){
      pop_diff <- utils.basic.stats(y)
      return(pop_diff$overall[c("Fst", "Fstp", "Dest", "Gst_H")])
    }
    x2 <- x[, indices]
    res_pop.diff_fun <- pop.diff_fun(x2)
    return(res_pop.diff_fun)
  }
  
  # DO THE JOB
  
  pops <- seppop(x)
  npops <- length(pops)
  pairs_pops <- t(combn(npops, 2))
  pairs_pops_names <- apply(pairs_pops, 1, function(y) {
    paste0(names(pops)[y[1]], "_vs_", names(pops)[y[2]])
  })
  
  ### pairwise
  if (npops > 2) {
    # observed value
    pairpop_res <- apply(pairs_pops, 1, function(y) {
      tpop <- rbind.dartR(pops[[y[1]]], pops[[y[2]]])
      res_tmp <-
        utils.basic.stats(tpop)$overall[c("Fst", "Fstp", "Dest", "Gst_H")]
      return(res_tmp)
    })
    
    if (nboots > 0) {
      # bootstrapping
      pairpop_boot <- apply(pairs_pops, 1, function(y) {
        tpop <- rbind.dartR(pops[[y[1]]], pops[[y[2]]])
        res_boots <- boot::boot(
          data = tpop,
          statistic = pop.diff,
          R = nboots,
          parallel = parallel,
          ncpus = ncpus
        )
        return(res_boots)
      })
      
      # confidence intervals
      # creating matrices to store CI
      res_CI <- replicate(length(pairpop_boot),
                          as.data.frame(matrix(nrow = 4, ncol = 2)),
                          simplify = FALSE)
      
      for (pop_n in 1:length(pairpop_boot)) {
        for (stat_n in 1:4) {
          res_CI_tmp <-     boot::boot.ci(
            boot.out = pairpop_boot[[pop_n]],
            conf = conf,
            type = CI.type,
            index = stat_n,
            t0 =  pairpop_res[stat_n, pop_n],
            t = pairpop_boot[[pop_n]]$t[, stat_n]
          )
          
          res_CI[[pop_n]][stat_n,] <-
            tail(as.vector(res_CI_tmp[[4]]), 2)
          
        }
      }
    }
    
  } else{
    tpop <- rbind.dartR(pops[[1]], pops[[2]])
    # observed values
    pairpop_res <-
      utils.basic.stats(tpop)$overall[c("Fst", "Fstp", "Dest", "Gst_H")]
    
    if (nboots > 0) {
      res_CI <- as.data.frame(matrix(nrow = 4, ncol = 2))
      
      # bootstrapping
      pairpop_boot <- boot::boot(
        data = tpop,
        statistic = pop.diff,
        R = nboots,
        parallel = parallel,
        ncpus = ncpus
      )
      
      # confidence intervals
      for (stat_n in 1:4) {
        
        res_CI_tmp <-
          boot::boot.ci(
            boot.out = pairpop_boot,
            conf = conf,
            type = CI.type,
            index = stat_n,
            t0 =  pairpop_res[stat_n],
            t = pairpop_boot$t[, stat_n]
          )
        
        res_CI[stat_n,] <-  tail(as.vector(res_CI_tmp[[4]]), 2)
        
      }
      
    }
  }
  
  if (npops > 2 & nboots > 0) {
    stat_pop <- asplit(pairpop_res, 2)
    stat_pop <- lapply(stat_pop, function(x) {
      as.data.frame(x)
    })
    CI <- Map(cbind, stat_pop, res_CI)
    CI <- lapply(CI, function(y) {
      colnames(y) <-  c("Value", "LCI", "HCI")
      return(y)
    })
    names(CI) <- pairs_pops_names
    
  }
  
  if (npops <= 2 & nboots > 0) {
    stat_pop <- pairpop_res
    CI <- cbind(stat_pop, res_CI)
    colnames(CI) <-  c("Value", "LCI", "HCI")
  }
  
  mat_pops <-
    rep(list(matrix(NA, nrow = npops, ncol = npops)), 4)
  
  if (npops > 2) {
    for (i in 1:length(mat_pops)) {
      mat_pops[[i]][lower.tri(mat_pops[[i]])] <- pairpop_res[i,]
      colnames(mat_pops[[i]]) <-
        rownames(mat_pops[[i]]) <- names(pops)
      mat_pops[[i]][upper.tri(mat_pops[[i]])] <-
        t(mat_pops[[i]])[rev(lower.tri(mat_pops[[i]]))]
    }
  } else{
    for (i in 1:length(mat_pops)) {
      mat_pops[[i]][lower.tri(mat_pops[[i]])] <- pairpop_res[i]
      colnames(mat_pops[[i]]) <-
        rownames(mat_pops[[i]]) <- names(pops)
      mat_pops[[i]][upper.tri(mat_pops[[i]])] <-
        t(mat_pops[[i]])[rev(lower.tri(mat_pops[[i]]))]
    }
  }
  
  if (npops > 2) {
    names(mat_pops) <- rownames(pairpop_res)
  } else{
    names(mat_pops) <- names(pairpop_res)
  }
  
  pairpop_res <- as.data.frame(pairpop_res)
  colnames(pairpop_res) <- pairs_pops_names
  
  # Printing outputs -----------
  if (verbose >= 2) {
    if (nboots > 0 & npops > 2) {
      print(list(
        Stat_matrices = mat_pops,
        Confident_Intervals = CI
      ))
    }
    
    if (nboots > 0 & npops <= 2) {
      print(list(
        Stat_tables = data.frame(pairpop_res),
        Confident_Intervals = CI
      ))
    }
    
    if (nboots == 0 & npops > 2) {
      print(list(
        Stat_matrices = mat_pops,
        Stat_tables = data.frame(Stat_tables = pairpop_res)
      ))
    }
    
    if (nboots == 0 & npops <= 2) {
      print(data.frame(Stat_tables = pairpop_res))
    }
  }
  
  # solution to print gplots https://stackoverflow.com/a/19191951
  create_heatmap <- function(...) {
    plot_heatmap <- function()
      gl.plot.heatmap(...)
  }
  
  p3 <- create_heatmap(
    mat_pops[[plot.stat]],
    palette.divergent = palette.divergent,
    cexRow = font.size,
    cexCol = font.size,
    na.color = "gray",
    symkey = FALSE,
    symbreaks = FALSE,
    verbose = verbose,
    ...
  )
  
  # PLOT THE RESULTS -----------------
  if (plot.display & npops > 2) {
    tryCatch(
      expr = {
        p3()
      },
      error = function(e) {
        cat(
          warn(
            "   Your plot was not shown in full because your 'Plots' pane
    is too small. Increase the size of the 'Plots' pane and run the
    function again. Alternatively, use the parameter 'plot.file' to
    save the plot to a file.\n"
          )
        )
      }
    )
    
  }
  
  if (plot.display & npops <= 2 & verbose >= 2) {
    cat(warn(
      "   No plot was displayed because only two populations were analysed.\n"
    ))
  }
  
  # Optionally save the plot ---------------------
  
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p3(),
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }
  
  # FINISH UP -------------------
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  # ----------------------
  
  # RETURN
  
  if (nboots > 0 & npops > 2) {
    return(list(
      Stat_matrices = mat_pops,
      Confident_Intervals = CI
    ))
  }
  
  if (nboots > 0 & npops <= 2) {
    return(list(Stat_tables = pairpop_res,
                Confident_Intervals = CI))
  }
  
  if (nboots == 0 & npops > 2) {
    return(list(Stat_matrices = mat_pops,
                data.frame(Stat_tables = pairpop_res)))
  }
  
  if (nboots == 0 & npops <= 2) {
    return(data.frame(Stat_tables = pairpop_res))
  }
  
}