#' @name gl.report.fstat
#' @title
#' Reports various statistics of genetic differentiation between
#' populations with confidence intervals
#' @family matched report
#' @description
#' This function calculates four genetic differentiation between populations
#' statistics (see the "Details" section for further information).
#'
#' \itemize{
#' \item \strong{Fst} - Measure of the degree of genetic differentiation of 
#' subpopulations (Nei, 1987).
#' \item \strong{Fstp} - Unbiased (i.e. corrected for sampling error, see 
#' explanation below) Fst (Nei, 1987).
#' \item \strong{Dest} - Jost's D (Jost, 2008).
#' \item \strong{Gst_H} - Gst standardized by the maximum level that it can obtain for
#' the observed amount of genetic variation (Hedrick 2005).
#' }
#' 
#' Sampling errors arise because allele frequencies in our samples differ from 
#' those in the subpopulations from which they were taken (Holsinger, 2012).
#'
#' Confidence Intervals are obtained by bootstrapping over loci.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param nboots Number of bootstrap replicates to obtain confidence intervals.
#' A whole number, either 0 (no bootstrap) or at least 2; with
#' CI.type = "bca", at least 200 [default 0].
#' @param conf The confidence level of the required interval, between 0 and 1
#' [default 0.95].
#' @param CI.type Method to estimate confidence intervals. One of
#' "norm", "basic", "perc" or "bca" [default "bca"].
#' @param ncpus Number of processes to be used in parallel operation. If ncpus
#' > 1 parallel operation is activated,see "Details" section [default 1].
#' @param plot.stat Statistic to plot. One of "Fst","Fstp","Dest" or "Gst_H"
#' [default "Fstp"].
#' @param plot.display If TRUE, a heatmap of the pairwise static chosen is
#'  displayed in the plot window [default TRUE].
#' @param palette.divergent A color palette function for the heatmap plot
#'  [default gl.colors("div")].
#' @param font.size Size of font for the labels of horizontal and vertical axes
#' of the heatmap [default 0.5].
#' @param plot.dir Directory in which to save files [default working directory].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
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
#'     The equations used to calculate the statistics are shown below.
#'
#'      \itemize{
#'      \item
#'     \emph{Ho} - Unbiased estimate of observed heterozygosity across 
#'     subpopulations (Nei, 1987, pp. 164, eq. 7.38) is calculated as:
#'
#'     \figure{Hoequation.jpg}
#'
#'     where \emph{Pkii} represents the proportion of homozygote \emph{ii} for 
#'     allele \emph{i} in individual \emph{k} and \emph{s} represents the number
#'      of subpopulations.
#'
#'     \item
#'     \emph{Hs} - Unbiased estimate of the expected heterozygosity under 
#'     Hardy-Weinberg equilibrium across subpopulations (Nei, 1987, pp. 164,
#'      eq. 7.39) is calculated as:
#'
#'     \figure{Hsequation.jpg}
#'     
#'     where \emph{ñ} is the harmonic mean of \emph{nk} (the number of 
#'     individuals in each subpopulation), \emph{pki} is the proportion 
#'     (sometimes misleadingly called frequency) of allele \emph{i} in 
#'     subpopulation \emph{k}. 
#'
#'     \item
#'     \emph{Ht} - Heterozygosity for the total population (Nei, 1987, pp. 164,
#'      eq. 7.40) is calculated as:
#'
#'     \figure{Htequation.jpg}
#'
#'     \item
#'     \emph{Dst} - The average allele frequency differentiation between 
#'     populations (Nei, 1987, pp. 163) is calculated as:
#'
#'     \figure{Dstequation.jpg}
#'
#'     \item
#'     \emph{Htp} - Unbiased estimate of Heterozygosity for the total population
#'     (Nei, 1987, pp. 165) is calculated as:
#'
#'     \figure{Htpequation.jpg}
#'
#'     \item
#'     \emph{Dstp} - Unbiased estimate of the average allele frequency
#'     differentiation between populations (Nei, 1987, pp. 165), where
#'     \emph{s} is the number of subpopulations that carry the locus, is
#'     calculated as:
#'
#'     \deqn{Dstp = Dst * s / (s - 1)}
#'
#'     \item
#'     \emph{Fst} - Measure of the extent of genetic differentiation 
#'     of subpopulations (Nei, 1987, pp. 162, eq. 7.34) is calculated as:
#'
#'     \figure{Fstequation.jpg}
#'
#'     \item
#'     \emph{Fstp} - Unbiased measure of the extent of genetic differentiation 
#'     of subpopulations (Nei, 1987, pp. 163, eq. 7.36) is calculated as:
#'
#'     \figure{Fstpequation.jpg}
#'
#'     \item
#'     \emph{Dest} - Jost's D (Jost, 2008, eq. 12) is calculated as:
#'
#'     \deqn{Dest = Dstp / (1 - Hs)}
#'
#'     \item
#'     \emph{Gst-max} - The maximum level that Gst can obtain for the observed 
#'     amount of genetic variation (Hedrick 2005, eq. 4a) is calculated as:
#'
#'     \figure{GstMaxequation.jpg}
#'
#'     \item
#'     \emph{Gst-H} - Gst standardized by the maximum level that it can obtain 
#'     for the observed amount of genetic variation (Hedrick 2005, eq. 4b) is 
#'     calculated as:
#'
#'     \figure{GstH.jpg}
#'
#'     }
#'
#'  \strong{Confidence Intervals}
#'
#' The uncertainty of a parameter, in this case the mean of the statistic, can
#' be summarised by a confidence interval (CI) which includes the true parameter
#' value with a specified probability (i.e. confidence level; the parameter
#' "conf" in this function).
#'
#' In this function, CI are obtained using Bootstrap which is an inference
#' method that samples with replacement the data and calculates the
#'  statistics every time.
#'
#'  The unit of resampling is the locus, which is the unit the statistics are
#'  averaged over. Each of the "nboots" replicates draws nLoc loci with
#'  replacement from the pair of populations being compared, and the four
#'  statistics are recalculated on that replicate by the same estimator that
#'  produces the reported value.
#'
#'  This function uses the function \link[boot]{boot} (package boot) to perform
#'  the bootstrap replicates and the function \link[boot]{boot.ci}
#'  (package boot) to perform the calculations for the CI.
#'
#'  The function has no seed parameter. To obtain the same intervals twice,
#'  set the global random number generator immediately before the call, for
#'  example \code{set.seed(1234)}.
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
#'    https://www.r-bloggers.com/2019/09/understanding-bootstrap-confidence-interval-output-from-the-r-boot-package/
#'
#'     Nice tutorials about the different types of CI can be found in:
#'
#'     https://www.datacamp.com/tutorial/bootstrap-r
#'
#'     and
#'
#'    https://www.r-bloggers.com/2019/09/understanding-bootstrap-confidence-interval-output-from-the-r-boot-package/
#'
#'      Efron and Tibshirani (1993, p. 162) and Davison and Hinkley
#'      (1997, p. 194) suggest that the number of bootstrap replicates should
#'      be between 1000 and 2000.
#'
#'  \strong{It is important} to note that unreliable confidence intervals will be
#'   obtained if too few number of bootstrap replicates are used.
#'   Therefore, the function \link[boot]{boot.ci} will throw warnings and errors
#'    if bootstrap replicates are too few. Consider increasing then number of
#'    bootstrap replicates to at least 200. With the default
#'    CI.type = "bca", fewer than 200 replicates is refused with an
#'    informative error rather than passed to \link[boot]{boot.ci}, which
#'    fails with "estimated adjustment 'a' is NA".
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
#' If a plot.file is given, the plot arising from this function is saved as an
#'  "RDS" binary file using the function \link[base]{saveRDS} (package base);
#'   can be reloaded with function \link[base]{readRDS} (package base). A file
#'   name must be specified for the plot to be saved.
#'
#'  If a plot directory (plot.dir) is specified, the gplot binary is saved to
#'  that directory; otherwise to the tempdir().
#'
#'  Your plot might not shown in full because your 'Plots' pane is too small
#'  (in RStudio).
#'  Increase the size of the 'Plots' pane before running the function.
#'  Alternatively, use the parameter 'plot.file' to save the plot to a file.
#'
#'  \strong{Parallelisation}
#'
#'  If the parameter ncpus > 1, parallelisation is enabled. In Windows, parallel
#'   computing employs a "socket" approach that starts new copies of R on each
#'    core. POSIX systems, on the other hand (Mac, Linux, Unix, and BSD),
#'    utilise a "forking" approach that replicates the whole current version of
#'     R and transfers it to a new core.
#'
#'     Opening and terminating R sessions in each core involves a significant
#'     amount of processing time, therefore parallelisation in Windows machines
#'    is only quicker than not using parallelisation when nboots > 1000-2000.
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
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
#' Holsinger, K. E. (2012). Lecture notes in population genetics.
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
#' @return The shape of the returned object depends on the number of
#' populations and on whether a bootstrap was requested. In every shape the
#' statistics are Fst, Fstp, Dest and Gst_H, and pairs of populations are
#' named "<pop1>_vs_<pop2>".
#' \itemize{
#' \item More than two populations, nboots > 0: a list of two elements.
#' "Stat_matrices" holds one nPop x nPop symmetric matrix per statistic;
#' "Confidence_Intervals" holds one data frame per pair of populations, with
#' one row per statistic and columns "Value", "LCI" (low confidence interval)
#' and "HCI" (high confidence interval).
#' \item Two populations, nboots > 0: a list of two elements. "Stat_tables"
#' is a data frame with one row per statistic and one column for the pair;
#' "Confidence_Intervals" is a single data frame with columns "Value", "LCI"
#' and "HCI".
#' \item More than two populations, nboots = 0: a list of two elements,
#' "Stat_matrices" as above and "Stat_tables", a data frame with one row per
#' statistic and one column per pair of populations.
#' \item Two populations, nboots = 0: a data frame with one row per statistic
#' and one column for the pair, not a list.
#' }
#' @export
#'
# ----------------------
# Function

gl.report.fstat <- function(x,
                            nboots = 0,
                            conf = 0.95,
                            CI.type = "bca",
                            ncpus = 1,
                            plot.stat = "Fstp",
                            plot.display = TRUE,
                            palette.divergent = gl.colors("div", verbose=0),
                            font.size = 0.5,
                            plot.dir = NULL,
                            plot.file = NULL,
                            verbose = NULL,
                            ...) {
  # PRELIMINARIES -- checking ----------------
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # verbose 0 is fully silent: no console output and no graphics [approved F4]
  if (verbose == 0) {
    plot.display <- FALSE
  }

  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   verbose = verbose)

  # CHECK DATATYPE
  # the arithmetic is SNP dosage specific throughout: heterozygotes are
  # scores of 1 and allele frequencies are the dosage mean halved, neither
  # of which means anything for presence/absence data [approved F3]
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING [approved F5, F7, F8]
  if (!is.numeric(nboots) || length(nboots) != 1 || is.na(nboots) ||
      nboots < 0 || nboots != round(nboots)) {
    stop(error(
      "Fatal Error: nboots must be a single non-negative whole number\n"))
  }

  if (nboots == 1) {
    stop(error(
      "Fatal Error: nboots = 1 gives a bootstrap distribution of a",
      "single value, from which no confidence interval can be",
      "calculated. Set nboots = 0 for point estimates alone, or at",
      "least 200 replicates for confidence intervals\n"))
  }

  if (length(CI.type) != 1 ||
      !CI.type %in% c("norm", "basic", "perc", "bca")) {
    stop(error(
      "Fatal Error: CI.type must be one of 'norm', 'basic', 'perc' or 'bca'\n"))
  }

  if (!is.numeric(conf) || length(conf) != 1 || is.na(conf) ||
      conf <= 0 || conf >= 1) {
    stop(error(
      "Fatal Error: conf must be a single number greater than 0",
      "and less than 1\n"))
  }

  if (length(plot.stat) != 1 ||
      !plot.stat %in% c("Fst", "Fstp", "Dest", "Gst_H")) {
    stop(error(
      "Fatal Error: plot.stat must be one of 'Fst', 'Fstp',",
      "'Dest' or 'Gst_H'\n"))
  }

  # boot.ci's bca interval needs enough replicates to estimate the
  # acceleration constant; below 200 it aborts with "estimated adjustment
  # 'a' is NA" from inside boot
  if (nboots > 0 && CI.type == "bca" && nboots < 200) {
    stop(error(
      "Fatal Error: CI.type = 'bca' requires at least 200 bootstrap",
      "replicates. Increase nboots, or choose CI.type = 'perc',",
      "'norm' or 'basic'\n"))
  }

  # keeping populations with more than 1 individuals
  pop_names <- popNames(x)[which(table(pop(x)) > 1)]

  if (length(pop_names) < 2) {
    stop(error(
      "Fatal Error: at least two populations of more than one",
      "individual each are required for pairwise comparisons\n"))
  }

  if (length(pop_names) < nPop(x)) {
    # dropping populations changes the result, so it is announced whenever
    # the user is listening at all [approved F10]
    if (verbose >= 1) {
      cat(warn(paste0(
        "  Keeping only populations with more than one individual. Dropped: ",
        paste(setdiff(popNames(x), pop_names), collapse = ", "),
        "\n"
      )))
    }
    x <- gl.keep.pop(x,
                     pop.list = pop_names,
                     verbose = verbose)
  }

  #converting to dartR object
  class(x) <- "dartR"
  
  # BOOTSTRAP STATISTIC [approved F1, F2] ----------
  # The unit of resampling is the locus, because that is the unit the four
  # statistics are averaged over. boot::boot draws its indices from the rows
  # of the data it is given, so it is handed a one-column data frame of
  # locus positions; the statistic then subsets the genotype matrix by
  # those positions. Previously boot was handed an individuals-by-loci
  # frame and the statistic applied the row indices to the columns, so only
  # the first nInd loci could ever enter a replicate.
  #
  # The replicate is rebuilt from the decoded genotype matrix with new()
  # rather than by subsetting the genlight, because adegenet's SNPbin "["
  # method drops NA on repeated indices and a draw with replacement always
  # produces repeated indices.
  #
  # The replicate statistics come from utils.basic.stats, the same function
  # that produces the reported point estimate, so the interval and the value
  # it brackets are the same estimator. The estimator was previously
  # duplicated inline here, and the copy had drifted.
  pop.diff <- function(loc.index,
                       indices,
                       gen.mat,
                       pops.info) {
    gl.boot <- new(
      "genlight",
      gen = gen.mat[, loc.index$loc[indices], drop = FALSE],
      ploidy = 2,
      pop = pops.info,
      parallel = FALSE
    )

    res <- utils.basic.stats(gl.boot)$overall[c("Fst", "Fstp", "Dest",
                                                "Gst_H")]

    return(res)

  }

  # setting parallel
  # if(ncpus>1){
  if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
    parallel <- "multicore"
  }
  ## if windows
  if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
    parallel <- "snow"
  }
  # }
  
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
          data = data.frame(loc = seq_len(nLoc(tpop))),
          statistic = pop.diff,
          gen.mat = as.matrix(tpop),
          pops.info = as.character(pop(tpop)),
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
          
          res_CI[[pop_n]][stat_n, ] <-
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
        data = data.frame(loc = seq_len(nLoc(tpop))),
        statistic = pop.diff,
        gen.mat = as.matrix(tpop),
        pops.info = as.character(pop(tpop)),
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
        
        res_CI[stat_n, ] <-  tail(as.vector(res_CI_tmp[[4]]), 2)
        
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
      mat_pops[[i]][lower.tri(mat_pops[[i]])] <- pairpop_res[i, ]
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
  # the results summary belongs at verbose 3; verbose 2 is a progress log
  # [approved F9]
  if (verbose >= 3) {
    if (nboots > 0 & npops > 2) {
      print(list(
        Stat_matrices = mat_pops,
        Confidence_Intervals = CI
      ))
    }

    if (nboots > 0 & npops <= 2) {
      print(list(
        Stat_tables = pairpop_res,
        Confidence_Intervals = CI
      ))
    }

    if (nboots == 0 & npops > 2) {
      print(list(
        Stat_matrices = mat_pops,
        Stat_tables = pairpop_res
      ))
    }

    if (nboots == 0 & npops <= 2) {
      print(pairpop_res)
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
        # report what actually failed; not every heatmap error is a pane
        # size problem [approved F7]
        cat(
          warn(
            paste0(
              "   The heatmap was not drawn. The error was: ",
              conditionMessage(e),
              "\n   If your 'Plots' pane is too small, increase its size and
    run the function again. Alternatively, use the parameter 'plot.file'
    to save the plot to a file.\n"
            )
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
      Confidence_Intervals = CI
    ))
  }
  
  if (nboots > 0 & npops <= 2) {
    return(list(Stat_tables = pairpop_res,
                Confidence_Intervals = CI))
  }
  
  # the pairwise table is named as a list element, not wrapped in a
  # data.frame() call that prefixes every column with "Stat_tables."
  # [approved F6, F11]
  if (nboots == 0 & npops > 2) {
    return(list(Stat_matrices = mat_pops,
                Stat_tables = pairpop_res))
  }

  if (nboots == 0 & npops <= 2) {
    return(pairpop_res)
  }
  
}