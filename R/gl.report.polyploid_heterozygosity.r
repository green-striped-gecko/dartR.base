#' @name gl.report.polyploid_heterozygosity
#' @title Reports observed, expected and unbiased heterozygosities and FIS
#' (inbreeding coefficient) by population or by individual from polyploid
#' (dosage) SNP data
#' @family unmatched report
#'
#' @description Calculates the observed (gametic), expected and unbiased
#' expected (i.e. corrected for sample size) heterozygosities and FIS
#' (inbreeding coefficient) for each population, or the observed
#' heterozygosity for each individual, in a genlight object coded as allele
#' dosages. For diploid data the results equal those of
#' \code{\link{gl.report.heterozygosity}}.

#' @param x Name of the genlight object containing SNP data coded as dosages
#' (0 to k copies of the alternative allele, where k is the ploidy), for
#' example read with \code{gl.read.vcf(mode = "dosage")} [required].
#' @param method Calculate heterozygosity by population (method='pop') or by
#' individual (method='ind') [default 'pop'].
#' @param n.invariant An estimate of the number of invariant sequence tags used
#' to adjust the heterozygosity rate [default 0].
#' @param subsample.pop Whether subsample populations to estimate observed 
#' heterozygosity (see Details) [default FALSE].
#' @param n.limit Minimum number of individuals that should have a population to 
#' perform subsampling to estimate heterozygosity [default 10].
#' @param nboots Number of bootstrap replicates to obtain confidence intervals
#' [default 0].
#' @param conf The confidence level of the required interval  [default 0.95].
#' @param CI.type Method to estimate confidence intervals. One of
#' "norm", "basic", "perc" or "bca" [default "bca"].
#' @param ncpus Number of processes to be used in parallel operation. If ncpus
#' > 1 parallel operation is activated, see "Details" section [default 1].
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors.pop A color palette for population plots or a list with
#' as many colors as there are populations in the dataset
#' [default gl.colors("dis")].
#' @param plot.colors.ind List of two color names for the borders and fill of
#' the plot by individual [default gl.colors(2)].
#' @param error.bar Statistic to be plotted as error bar either "SD" (standard 
#' deviation) or "SE" (standard error) or "CI" (confidence intervals)
#'  [default "SD"].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()].
#' @param plot.file Name for the RDS binary file to save (base name only, 
#' exclude extension) [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#'
#' @details
#' Observed heterozygosity is the gametic heterozygosity: for an individual
#' of ploidy k carrying d copies of the alternative allele, the probability
#' that two allele copies drawn without replacement differ,
#' d(k - d) / choose(k, 2) (Moody et al. 1993). For diploids this is 1 for a
#' heterozygote and 0 for a homozygote. Observed heterozygosity for a
#' population averages it over the individuals scored at each locus and then
#' over loci. The calculations take into account missing values.
#'
#' Expected heterozygosity for a population takes the expected proportion of
#' heterozygotes, that is, expected under Hardy-Weinberg equilibrium, for each
#' locus, then averages this across the loci for an average estimate for the
#' population.
#'
#' The unbiased expected heterozygosity is calculated using the correction for 
#' sample size following equation 2 from Nei 1978.
#' 
#' Accuracy of all heterozygosity estimates is affected by small sample sizes,
#' and so is their comparison between populations or repeated analysis. Expected
#' heterozygosities are less affected because their calculations are based on 
#' allele frequencies while observed heterozygosities are strongly susceptible 
#' to sampling effects when the sample size is small.  
#'
#' Observed heterozygosity for individuals is the gametic heterozygosity
#' averaged over the loci scored for that individual. The output also gives
#' the proportions of loci homozygous for the reference allele (dosage 0) and
#' for the alternative allele (dosage k).
#'
#' Finally, the loci that are invariant across all individuals in the dataset
#' (that is, across populations), is typically unknown. This can render
#' estimates of heterozygosity analysis specific, and so it is not valid to
#' compare such estimates across species or even across different analyses 
#' (see Schmidt et al 2021). This is a similar problem faced by microsatellites. 
#' If you have an estimate of the
#' number of invariant sequence tags (loci) in your data, such as provided by
#' \code{\link{gl.report.secondaries}}, you can specify it with the n.invariant
#' parameter to standardize your estimates of heterozygosity. This is called
#' autosomal heterozygosities by Schmidt et al (2021).
#'
#' \strong{NOTE}: It is important to realise that estimation of adjusted (autosomal)
#' heterozygosity requires that secondaries not to be removed.
#'
#' Heterozygosities and FIS (inbreeding coefficient) are calculated by locus
#' within each population using the following equations, and then averaged across 
#' all loci:
#' \itemize{
#' \item Observed heterozygosity (Ho) = mean over individuals of
#' d(k - d) / choose(k, 2), where d is the dosage and k the ploidy.
#' \item Observed heterozygosity adjusted (Ho.adj) <- Ho * n_Loc /
#'  (n_Loc + n.invariant),
#' where n_Loc is the number of loci that do not have all missing data  and
#' n.invariant is an estimate of the number of invariant loci to adjust
#' heterozygosity.
#' \item Expected heterozygosity (He) = 1 - (p^2 + q^2),
#' where q, the frequency of the alternative allele, is the sum of dosages
#' divided by the number of sampled allele copies (the sum of ploidies of the
#' individuals scored at that locus), and p = 1 - q.
#' \item Expected heterozygosity adjusted (He.adj) = He * n_Loc /
#' (n_Loc + n.invariant)
#' \item Unbiased expected heterozygosity (uHe) = He * N / (N - 1), where N
#' is the number of sampled allele copies at that locus (2 * n_Ind for
#' diploids)
#' \item Inbreeding coefficient (FIS) = 1 - Ho / uHe
#' }

#'\strong{ Function's output }
#'
#' Output for method='pop' is an ordered barchart of observed heterozygosity,
#' unbiased expected heterozygosity and FIS (Inbreeding coefficient) across 
#' populations together with a table of mean observed and expected 
#' heterozygosities and FIS by population and their respective standard 
#' deviations (SD).

#' In the output, it is also reported by population: the number of loci used to
#'  estimate heterozygosity (n.Loc), the number of polymorphic loci (polyLoc),
#'  the number of monomorphic loci (monoLoc) and loci with all missing data
#'   (all_NALoc).
#'   
#' Output for method='ind' is a histogram and a boxplot of heterozygosity across
#' individuals.
#' 
#' If a plot.file is given, the ggplot arising from this function is saved as an 
#' "RDS" binary file using saveRDS(); can be reloaded with readRDS(). A file 
#' name must be specified for the plot to be saved.
#' 
#' If a plot directory (plot.dir) is specified, the ggplot binary is saved to 
#' that directory; otherwise to the tempdir(). Nothing is saved when
#' plot.display = FALSE (or verbose = 0), as no plot is built.
#'  
#'  Examples of other themes that can be used can be consulted in: 
#'  \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'  
#'  \strong{Subsampling populations}
#'  
#' Subsampling applies to method = 'pop' only; it is ignored, with a warning,
#' for method = 'ind'. To test the effect of five population sample sizes (n = 10, 5, 4, 3, 2) on 
#' observed heterozygosity estimates, the function subsamples individuals,
#'  without replacement. The subsampling is repeated 10 times for each sample
#'   size n. This approach is an implementation of Schmidt et al (2021). 
#'  
#'   \strong{Error bars}
#'  
#'  The best method for presenting or assessing genetic statistics depends on 
#'  the type of data you have and the specific questions you're trying to 
#'  answer. Here's a brief overview of when you might use each method:
#'  
#'   \strong{1. Confidence Intervals ("CI"):}
#'   
#'  - Usage: Often used to convey the precision of an estimate.
#'  
#'  - Advantage: Confidence intervals give a range in which the true parameter 
#'  (like a population mean) is likely to fall, given the data and a specified 
#'  probability (like 95\%).
#'  
#'  - In Context: For genetic statistics, if you're estimating a parameter,
#'   a 95\% CI gives you a range in which you're 95\% confident the true parameter
#'    lies.
#'  
#'   \strong{2. Standard Deviation ("SD"):}
#'   
#'  - Usage: Describes the amount of variation from the average in a set of data.
#'  
#'  - Advantage: Allows for an understanding of the spread of individual data
#'   points around the mean.
#'   
#'  - In Context: If you're looking at the distribution of a quantitative trait 
#'  (like height) in a population with a particular genotype, the SD can 
#'  describe how much individual heights vary around the average height.
#'  
#'   \strong{3. Standard Error ("SE"):}
#'   
#'  - Usage: Describes the precision of the sample mean as an estimate of the 
#'  population mean.
#'  
#'  - Advantage: Smaller than the SD in large samples; it takes into account 
#'  both the SD and the sample size. 
#'  
#'  - In Context: If you want to know how accurately your sample mean represents
#'   the population mean, you'd look at the SE.
#'   
#'    \strong{Recommendation:}
#'    
#'   - If you're trying to convey the precision of an estimate, confidence 
#'   intervals are very useful.
#'   
#'   - For understanding variability within a sample, standard deviation is key.
#'   
#'   - To see how well a sample mean might estimate a population mean, consider 
#'   the standard error.
#'   
#'  \strong{Confidence Intervals}
#'
#' The uncertainty of a parameter, in this case the mean of the statistic, can
#' be summarised by a confidence interval (CI) which includes the true parameter
#' value with a specified probability (i.e. confidence level; the parameter
#' "conf" in this function).
#'
#' In this function, CI are obtained using Bootstrap which is an inference
#' method that samples loci with replacement and calculates the statistics
#'  every time.
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
#'    https://www.r-bloggers.com/2019/09/understanding-bootstrap-confidence-interval-output-from-the-r-boot-package
#'
#'     Nice tutorials about the different types of CI can be found in:
#'
#'     https://www.datacamp.com/tutorial/bootstrap-r
#'
#'     and
#'
#'    https://www.r-bloggers.com/2019/09/understanding-bootstrap-confidence-interval-output-from-the-r-boot-package
#'
#'      Efron and Tibshirani (1993, p. 162) and Davison and Hinkley
#'      (1997, p. 194) suggest that the number of bootstrap replicates should
#'      be between 1000 and 2000.
#'
#'  \strong{It is important} to note that unreliable confidence intervals will be
#'   obtained if too few number of bootstrap replicates are used.
#'   Therefore, the function \link[boot]{boot.ci} will throw warnings and errors
#'    if bootstrap replicates are too few. Consider increasing the number of
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
#' @author Author(s): Ching Ching Lau. Custodian: Ching Ching Lau (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' 
#'
#' @references
#' \itemize{
#' \item Moody, M. E., Mueller, L. D., & Soltis, D. E. (1993). 
#' Genetic variation and random drift in autotetraploid populations. Genetics, 134(2), 649-657.
#' \item Nei, M. (1978). Estimation of average heterozygosity and genetic
#' distance from a small number of individuals. Genetics, 89(3), 583-590.
#' \item Schmidt, T. L., Jasper, M., Weeks, A. R., & Hoffmann, A. A. (2021).
#' Unbiased population heterozygosity estimates from genome-wide sequence
#' data. Methods in Ecology and Evolution, 12(10), 1888-1898.
#'   }
#'
#' @examples
#' # simulated autotetraploid: 20 individuals in two populations, 50 loci,
#' # dosages 0-4 copies of the alternative allele
#' set.seed(1)
#' q <- runif(50, 0.1, 0.9)
#' dos <- sapply(q, function(qq) rbinom(20, 4, qq))
#' rownames(dos) <- paste0("ind", 1:20)
#' colnames(dos) <- paste0("loc", 1:50)
#' tetra <- new("genlight", dos, ploidy = rep(4, 20))
#' pop(tetra) <- rep(c("A", "B"), each = 10)
#' tetra <- gl.compliance.check(tetra, verbose = 0)
#' res <- gl.report.polyploid_heterozygosity(tetra)
#' res_ind <- gl.report.polyploid_heterozygosity(tetra, method = "ind")

#' @seealso \code{\link{gl.filter.heterozygosity}}

#' @export
#' @return For method = 'pop', a dataframe containing population labels,
#' heterozygosities, FIS, their standard deviations, standard errors,
#' confidence intervals (if nboots > 0) and sample sizes. For method = 'ind',
#' a dataframe with, per individual, Ho, the proportions of homozygous
#' reference and alternative loci, and the number of loci scored. With
#' subsample.pop = TRUE (method = 'pop'), a named list:
#' \code{subsample}, the subsampling results, and \code{results}, the
#' dataframe above.

gl.report.polyploid_heterozygosity <- function(x,
                                     method = "pop",
                                     n.invariant = 0,
                                     subsample.pop = FALSE,
                                     n.limit = 10,
                                     nboots = 0,
                                     conf = 0.95,
                                     CI.type = "bca",
                                     ncpus = 1,
                                     plot.display = TRUE,
                                     plot.theme = theme_dartR(),
                                     plot.colors.pop = gl.colors("dis"),
                                     plot.colors.ind = gl.colors(2),
                                     plot.file = NULL,
                                     plot.dir = NULL,
                                     error.bar = "SD",
                                     verbose = NULL) {
  
  # Utils functions are in utils.het.report.r
  
  # setting parallel
  # if (ncpus > 1) {
    if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
      parallel <- "multicore"
    }
    ## if windows
    if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
      parallel <- "snow"
    }
  # }
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  if(verbose==0){plot.display <- FALSE}

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jackson",
                   verbose = verbose)

  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir,verbose=0)

  # CHECK DATATYPE
  datatype <-
    utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # Gametic heterozygosity needs at least two allele copies per individual
  if (any(ploidy(x) < 2)) {
    stop(error(
      "  Fatal Error: all individuals must have ploidy >= 2 (dosage data)\n"
    ))
  }
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  if (!(method == "pop" | method == "ind")) {
    if (verbose >= 1) {
      cat(
        warn(
          "Warning: Method must either be by population or by individual,
                set to method='pop'\n"
        )
      )
    }
    method <- "pop"
  }

  if (n.invariant < 0) {
    if (verbose >= 1) {
      cat(warn(
        "Warning: Number of invariant loci must be non-negative, set to
            zero\n"
      ))
    }
    n.invariant <- 0
    if (verbose == 5) {
      cat(report(
        "  No. of invariant loci can be esimated using
                    gl.report.secondaries\n"
      ))
    }
  }

  if (any(grepl(x@other$history, pattern = "gl.filter.secondaries") == TRUE) &
      n.invariant > 0) {
    if (verbose >= 1) {
      cat(
        warn(
          "  Warning: Estimation of adjusted heterozygosity requires that
                secondaries not to be removed. A gl.filter.secondaries call was
                found in the history. This may cause the results to be
                incorrect\n"
        )
      )
    }
  }

  if( nboots == 0 & error.bar == "CI" ){
    stop(error(
      "  Number of bootstraps ('nboots' parameter) must be > 0 to calculate confidence intervals\n"
    ))
  }

  if (subsample.pop == TRUE && method == "ind") {
    if (verbose >= 1) {
      cat(warn(
        "  Warning: subsample.pop applies to method='pop' only -- ignored\n"
      ))
    }
    subsample.pop <- FALSE
  }

  # DO THE JOB
  
  ########### FOR METHOD BASED ON POPULATIONS
  
  if (method == "pop") {
    # Set a population if none is specified (such as if the genlight object
    # has been generated manually)
    if (is.null(pop(x)) |
        is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2) {
        cat(
          warn(
            "  No population assignments detected,
                             individuals assigned to a single population
                        labelled 'pop1'\n"
          )
        )
      }
      pop(x) <- array("pop1", dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }
    
    # Subsampling populations as in Schmidt (2021)
    if(subsample.pop==TRUE){
    res_sub <- utils.subsample.pop(x,n.limit = n.limit)
    }
    
    # Split the genlight object into a list of populations
    sgl <- seppop(x)
    
    # One genotype matrix and one ploidy vector per population
    sgl_m <- lapply(sgl, as.matrix)
    sgl_k <- lapply(sgl, function(y) as.numeric(ploidy(y)))

    # Calculate the number of individuals
    n_ind <- sapply(sgl, ind.count)

    # Number of polymorphic, monomorphic and all-NA loci by population. A
    # locus is monomorphic when every sampled allele copy is the same allele
    # (alternative-allele share 0 or 1).
    loc_counts <- mapply(function(m, k) {
      all_na <- colSums(!is.na(m)) == 0
      q <- colSums(m, na.rm = TRUE) / colSums((!is.na(m)) * k)
      mono <- !all_na & (q == 0 | q == 1)
      c(poly = sum(!all_na & !mono), mono = sum(mono), all_na = sum(all_na))
    }, sgl_m, sgl_k)
    poly_loc <- unname(loc_counts["poly", ])
    mono_loc <- unname(loc_counts["mono", ])
    all_na_loc <- unname(loc_counts["all_na", ])
    # Calculate the number of loci that are not all NAs CP ###
    n_loc <- poly_loc + mono_loc


    #### Calculate heterozygosities for each population ####
    if (verbose >= 2) {
      cat(
        report(
          "  Calculating Heterozygosities, averaged across
                    loci, for each population\n"
        )
      )
    }
    
    all.het <- mapply(function(m, k) pop.het_fun(m,
                                                  n.invariant = n.invariant,
                                                  aHet = n.invariant > 0,
                                                  bootstrap = FALSE,
                                                  ploidy = k),
                      sgl_m, sgl_k, SIMPLIFY = FALSE)
    
    Ho.loc <- sapply(all.het, function(x) x[["byloc"]]["Ho.loc"])
    
    He.loc <- sapply(all.het, function(x) x[["byloc"]]["He.loc"])
      
    Ho <- sapply(all.het, function(x) x[["means"]]["Ho"])
    HoSD <- compute.variability(all.het, what.st = "sd", what.het = "Ho.loc")
    HoSE <- compute.variability(all.het, what.st = "std.error", what.het = "Ho.loc")
      
    Hexp <- sapply(all.het, function(x) x[["means"]]["He"])
    HexpSD <- compute.variability(all.het, what.st = "sd", what.het = "He.loc")
    HexpSE <- compute.variability(all.het, what.st = "std.error", what.het = "He.loc")
    
    uHexp <- sapply(all.het, function(x) x[["means"]]["uHe"])
    uHexpSD <- compute.variability(all.het, what.st = "sd", what.het = "uHe.loc")
    uHexpSE <- compute.variability(all.het, what.st = "std.error", what.het = "uHe.loc")
    
    FIS <- sapply(all.het, function(x) x[["means"]]["FIS"])
    FISSD <- compute.variability(all.het, what.st = "sd", what.het = "FIS.loc")
    FISSE <- compute.variability(all.het, what.st = "std.error", what.het = "FIS.loc")
    
    
    if (n.invariant > 0) {
      # Apply correction CP ###
      Ho.adj <- sapply(all.het, function(x) x[["means"]]["Ho.adj"])
      # Manually compute SD for Ho.adj sum of the square of differences from
      # the mean for polymorphic sites plus sum of the square of differences
      #(which is the Ho.adj because Ho=0) from the mean for invariant sites
      Ho.adjSD <-
        sqrt((
          mapply(function(x, Mean)
            sum((x - Mean) ^ 2, na.rm = TRUE), Ho.loc, Mean = Ho.adj) +
            n.invariant * Ho.adj ^ 2) / (n_loc + n.invariant - 1))
      
      Ho.adjSE <-  Ho.adjSD / sqrt(n_loc+n.invariant)
      
      Hexp.adj <- sapply(all.het, function(x) x[["means"]]["Hexp.adj"])
      Hexp.adjSD <- 
        sqrt((
          mapply(function(x, Mean)
            sum((x - Mean) ^ 2, na.rm = TRUE), He.loc, Mean = Hexp.adj) + 
            n.invariant * Hexp.adj ^ 2) /(n_loc + n.invariant - 1))
      Hexp.adjSE <- Hexp.adjSD / sqrt(n_loc+n.invariant)
    }
    
    # Prep for results
    df.base <-
      data.frame(
        pop = popNames(x),
        n.Ind = round(n_ind, 6),
        n.Loc = n_loc,
        n.Loc.adj = n_loc / (n_loc + n.invariant),
        polyLoc = poly_loc ,
        monoLoc = mono_loc ,
        all_NALoc = all_na_loc)
    
    if (n.invariant == 0) {
      df.params <- data.frame(
        Ho = round(as.numeric(Ho),6),
        HoSD = round(HoSD,6),
        HoSE = round(HoSE, 6),
            
        He = round(Hexp, 6),
        HeSD = round(HexpSD, 6),
        HeSE = round(HexpSE,6),
       
        uHe = round(uHexp, 6),
        uHeSD = round(uHexpSD, 6),
        uHeSE = round(uHexpSE ,6),
       
        FIS = round(FIS,6),
        FISSD = round(FISSD,6),
        FISSE = round(FISSE , 6)
    )
    } else {
      df.params <- data.frame(
        Ho.adj = round(as.numeric(Ho.adj),6),
        Ho.adjSD = round(Ho.adjSD,6),
        Ho.adjSE = round(Ho.adjSE, 6),
        
        He.adj = round(Hexp.adj, 6),
        He.adjSD = round(Hexp.adjSD, 6),
        He.adjSE = round(Hexp.adjSE, 6)
        )
    }
    
    df <- cbind(df.base, df.params)
    npops <- nPop(x)
    
    #### bootstrapping ####
    if (nboots > 0) {
      # Loci are resampled: boot() resamples rows, so the matrix is passed
      # as loci x individuals and pop.het() transposes it back
      pop_boot <- mapply(function(m, k) {
        boot::boot(
          data = t(m),
          statistic = pop.het,
          n.invariant = n.invariant,
          aHet = n.invariant > 0,
          boot_method = "loc",
          ploidy = k,
          R = nboots,
          parallel = parallel,
          ncpus = ncpus
        )
      }, sgl_m, sgl_k, SIMPLIFY = FALSE)
      
      # confidence intervals
      
      # creating matrices to store CI
      nparams <- ncol(pop_boot[[1]]$t)
      res_CI <- replicate(npops,
                          as.data.frame(matrix(nrow = nparams, ncol = 2)),
                          simplify = FALSE)
      
      if(nparams == 4) {
        pop_res <- rbind(Ho, Hexp, uHexp, FIS)
      } else {
        pop_res <- rbind(Ho.adj, Hexp.adj)
      }
      
      no_CI <- character(0)
      for (pop_n in 1:length(sgl)) {
        for (stat_n in seq_len(nparams)) {
          # boot.ci cannot produce an interval when every replicate is
          # identical (e.g. a single-individual population resampled by
          # individual; it returns NULL) or when replicates are NA (e.g.
          # FIS where uHe is 0; it errors). Record NA limits and report
          # the population below instead of aborting the whole run.
          # (capture.output: boot.ci print()s its constant-t message,
          # which would break silence at verbose = 0)
          res_CI_tmp <- tryCatch({
            utils::capture.output(
              ci <- suppressWarnings(boot::boot.ci(
                boot.out = pop_boot[[pop_n]],
                conf = conf,
                type = CI.type,
                index = stat_n
              ))
            )
            ci
          }, error = function(e) NULL)
          if (is.null(res_CI_tmp) || length(res_CI_tmp) < 4) {
            no_CI <- union(no_CI, names(sgl)[pop_n])
            res_CI[[pop_n]][stat_n,] <- c(NA_real_, NA_real_)
          } else {
            res_CI[[pop_n]][stat_n,] <-
              tail(as.vector(res_CI_tmp[[4]]), 2)
          }
          
        }
      }
      if (length(no_CI) > 0 && verbose >= 1) {
        cat(warn(
          "  Warning: no confidence interval could be computed for",
          length(no_CI), "population(s); limits set to NA:",
          paste(no_CI, collapse = ", "), "\n"
        ))
      }
      
      # Build df with CI
      params_names <- names(df.params)[grep(pattern = "SD|SE", x = names(df.params), invert = TRUE)]
      lCI_df <- lapply(seq_len(nparams), function(rn, lres=res_CI, nms=params_names) {
        df <- do.call(rbind, lapply(lres, "[", rn,))
        names(df) <- paste0(params_names[rn], c("LCI", "HCI"))
        return(df)
      })
      df.CI <- do.call(cbind, lCI_df)
      df <- cbind(df, df.CI)
      
      # set up name order for df
      sfx <- c("", "SD", "SE", "LCI", "HCI")
      nms.order <- vector("character", length = nparams * length(sfx))
      for(rn in seq_len(nparams)) {
        nms.order[seq_len(length(sfx)) + length(sfx) * (rn - 1)] <- 
          paste0(params_names[rn], sfx)
      }
      # cbind CI to existing results
      df <- df[, c(names(df.base), nms.order)]
    }
    
    if (plot.display) {
      res.mean <- subsample <- error_L <- error_H <- value <- color <- variable <- He.adj <- res_SE <- NULL
      
      pop_order <- unique(as.character(pop(x))) 
      # printing plots and reports assigning colors to populations
      if (is(plot.colors.pop, "function")) {
        colors_pops <- plot.colors.pop(length(levels(pop(x))))
      }
      
      if (!is(plot.colors.pop, "function")) {
        colors_pops <- plot.colors.pop
      }
      colors_pops <- setNames(colors_pops, pop_order)
      
      if (n.invariant == 0) {
        
        pop_list_plot <- df
        pop_list_plot$pop <- as.factor(pop_list_plot$pop)
        pop_list_plot$color <- colors_pops
        
        pop_list_plot_stat <- pop_list_plot[,c("Ho", "uHe", "FIS", "n.Ind",  "pop",  "color")]
        pop_list_plot_stat <- reshape2::melt(pop_list_plot_stat, id = c("pop", "color", "n.Ind"))
        
        if(error.bar=="SD"){
          pop_list_plot_error <- pop_list_plot[,c("HoSD", "uHeSD", "FISSD","pop")]
          pop_list_plot_error <- reshape2::melt(pop_list_plot_error,id = c("pop"))
          colnames(pop_list_plot_error) <- c("pop","variable","error")
          pop_list_plot_error <- pop_list_plot_error[,c("pop","error")]
          pop_list_plot_stat <- cbind(pop_list_plot_stat,error=pop_list_plot_error$error)
        }
        
        if(error.bar=="SE"){
          pop_list_plot_error <- pop_list_plot[,c("HoSE","uHeSE","FISSE","pop")]
          pop_list_plot_error <- reshape2::melt(pop_list_plot_error,id = c("pop"))
          colnames(pop_list_plot_error) <- c("pop","variable","error")
          pop_list_plot_error <- pop_list_plot_error[,c("pop","error")]
          pop_list_plot_stat <- cbind(pop_list_plot_stat,error=pop_list_plot_error$error)
          }
        
        if(error.bar=="CI"){
          pop_list_plot_error_L <- pop_list_plot[,c("HoLCI","uHeLCI","FISLCI","pop")]
          pop_list_plot_error_H <- pop_list_plot[,c("HoHCI","uHeHCI","FISHCI","pop")]
          pop_list_plot_error_L <- reshape2::melt(pop_list_plot_error_L,id = c("pop"))
          pop_list_plot_error_H <- reshape2::melt(pop_list_plot_error_H,id = c("pop"))
          colnames(pop_list_plot_error_L) <- c("pop","variable_L","error_L")
          colnames(pop_list_plot_error_H) <- c("pop","variable_H","error_H")
          pop_list_plot_error <- cbind(pop_list_plot_error_L,pop_list_plot_error_H)
          pop_list_plot_stat <- cbind(pop_list_plot_stat,
                                      error_L = pop_list_plot_error$error_L,
                                      error_H = pop_list_plot_error$error_H)
          
          }
        
        pop_list_plot_stat$pop <- factor(pop_list_plot_stat$pop, levels = pop_order)
        
        lab_df <- pop_list_plot_stat[!duplicated(pop_list_plot_stat$pop),
                                     c("pop","n.Ind")]
        labels_named <- setNames(paste(lab_df$pop, round(lab_df$n.Ind, 0), sep = " | "),
                                 lab_df$pop)
        p3 <-
          ggplot(data = pop_list_plot_stat, aes(x = pop, 
                                                y = value,
                                                fill = pop)) +
          geom_bar(stat = "identity", 
                   color = "black", 
                   position = position_dodge())+ 
          facet_wrap(~variable, nrow=1) +
          scale_fill_manual(values = colors_pops,
                             breaks = pop_order,
                             limits = pop_order) +
          scale_x_discrete(limits = pop_order, labels = labels_named) +
          # scale_x_discrete(labels = paste(pop_list_plot_stat$pop,
          #                                 round(pop_list_plot_stat$n.Ind,
          #                                       0),
          #                                 sep = " | ")) +
          plot.theme +
          theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(
              angle = 90,
              hjust = 1,
              face = "bold",
              size = 12
            ),
            axis.title.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none"
          ) 
        
        if(error.bar=="SD"){
          p3 <- p3 + 
            geom_errorbar(aes(ymin = value, 
                              ymax = value + error), 
                          width=0.5)+
            ggtitle(label = "Heterozygosities and FIS by Population",
                    subtitle = "Error bars show Standard Deviation")
        }
        
        if(error.bar=="SE"){
          p3 <- p3 + 
            geom_errorbar(aes(ymin = value - error, 
                              ymax = value + error), 
                          width=0.5) +
            ggtitle(label = "Heterozygosities and FIS by Population",
                    subtitle = "Error bars show Standard Error")
        }
        
        if(error.bar=="CI"){
          p3 <- p3 + 
            geom_errorbar(aes(ymin = error_L, 
                              ymax = error_H), 
                          width=0.5)+
            ggtitle(label = "Heterozygosities and FIS by Population",
                    subtitle = "Error bars show Confidence Intervals")
        }
        
        # Subsampling populations as in Schmidt (2021)
        if(subsample.pop==TRUE){
          
          res_sub_plot <- res_sub
          res_sub_plot$pop <- as.factor(res_sub_plot$pop)
          res_sub_plot$subsample <- as.factor(res_sub_plot$subsample )
          # key colours by population name; populations below n.limit are
          # skipped by utils.subsample.pop, so a positional rep() desyncs
          res_sub_plot$color <- colors_pops[as.character(res_sub_plot$pop)]
          
          res_sub_plot_2 <- reshape2::melt(res_sub_plot, id = c("pop", "color", "subsample","res_SE"))
          
          p4 <- ggplot(res_sub_plot_2,aes(x=subsample,y=value))+
            
            geom_bar(stat = "identity", 
                     color = "black", 
                     fill = res_sub_plot_2$color,
                     position = position_dodge())+ 
            facet_wrap(~variable, nrow=1) +
            facet_wrap(~pop)+
            geom_errorbar(aes(ymin = value, 
                              ymax = value + res_SE), 
                          width=0.5)+
            ggtitle(label = "Observed heterozygosity in populations subsamples",
                    subtitle = "Error bars show Standard Error") +
            plot.theme +
            theme(
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(
                face = "bold",
                size = 12
              ),
              axis.title.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none"
            )
          
          p3 <-  p3 / p4
          
        }
        
      } else {

        df.ordered <- df
        df.ordered$color <- colors_pops
        df.ordered <- df.ordered[order(df.ordered$Ho.adj), ]
        df.ordered$pop <- factor(df.ordered$pop, levels = df.ordered$pop)
        p1 <-
          ggplot(df.ordered, aes(
            x = pop,
            y = Ho.adj,
            fill = pop
          )) + geom_bar(position = "dodge",
                        stat = "identity",
                        color = "black") +
          scale_fill_manual(values = df.ordered$color) +
          scale_x_discrete(labels = paste(df.ordered$pop,
                                          round(df.ordered$n.Ind, 0),
                                          sep = " | ")) + plot.theme + theme(
                                            axis.ticks.x = element_blank(),
                                            axis.text.x = element_blank(),
                                            axis.title.x = element_blank(),
                                            axis.ticks.y = element_blank(),
                                            axis.title.y = element_blank(),
                                            legend.position = "none"
                                          ) +
          labs(fill = "Population") +
          ggtitle("Adjusted Observed Heterozygosity by Population")
        
        p2 <-
          ggplot(df.ordered, aes(
            x = pop,
            y = He.adj,
            fill = pop
          )) + geom_bar(position = "dodge",
                        stat = "identity",
                        color = "black") +
          scale_fill_manual(values = df.ordered$color) +
          scale_x_discrete(labels = paste(df.ordered$pop,
                                          round(df.ordered$n.Ind, 0),
                                          sep = " | ")) + plot.theme + theme(
                                            axis.ticks.x = element_blank(),
                                            axis.text.x = element_text(
                                              angle = 90,
                                              hjust = 1,
                                              face = "bold",
                                              size = 12
                                            ),
                                            axis.title.x = element_blank(),
                                            axis.ticks.y = element_blank(),
                                            axis.title.y = element_blank(),
                                            legend.position = "none"
                                          ) +
          labs(fill = "Population") +
          ggtitle("Adjusted Expected Heterozygosity by Population")
        
        p3 <- (p1 / p2)
      }
    }
    
    # OUTPUT REPORT
    if (verbose >= 3) {
      cat("  Reporting Heterozygosity by Population\n")
      cat("\n  No. of loci =", nLoc(x), "\n")
      cat("  No. of individuals =", nInd(x), "\n")
      cat("  No. of populations =", nPop(x), "\n")
      if (n.invariant == 0) {
      cat("    Minimum Observed Heterozygosity: ", round(min(df$Ho, na.rm = TRUE), 6), "\n")
      cat("    Maximum Observed Heterozygosity: ", round(max(df$Ho, na.rm = TRUE), 6), "\n")
      cat("    Average Observed Heterozygosity: ", round(mean(df$Ho, na.rm = TRUE), 6), "\n\n")
      
      cat("    Minimum Unbiased Expected Heterozygosity: ",
          round(min(df$uHe, na.rm = TRUE), 6), "\n")
      cat("    Maximum Unbiased Expected Heterozygosity: ",
          round(max(df$uHe, na.rm = TRUE), 6), "\n")
      cat("    Average Unbiased Expected Heterozygosity: ",
          round(mean(df$uHe, na.rm = TRUE), 6), "\n")
      cat("  Heterozygosity estimates not corrected for uncalled invariant loci\n")
      
      } else {
        
        cat("    Minimum Observed adjusted Heterozygosity: ", 
            round(min(df$Ho.adj, na.rm = TRUE), 6), "\n")
        cat("    Maximum Observed adjusted Heterozygosity: ", 
            round(max(df$Ho.adj, na.rm = TRUE), 6), "\n")
        cat("    Average Observed adjusted Heterozygosity: ",
            round(mean(df$Ho.adj, na.rm = TRUE), 6), "\n\n")
        cat("    Minimum Unbiased adjusted Expected Heterozygosity: ", 
            round(min(df$He.adj, na.rm = TRUE), 6), "\n")
        cat("    Maximum Unbiased adjusted Expected Heterozygosity: ", 
            round(max(df$He.adj, na.rm = TRUE), 6), "\n")
        cat("    Average Unbiased adjusted Expected Heterozygosity: ",
            round(mean(df$He.adj, na.rm = TRUE), 6), "\n\n")
        cat(
          "  Average correction factor for invariant loci =",
          mean(n_loc / (n_loc + n.invariant), na.rm = TRUE),
          "\n"
        )
      }
    }
      
    # PRINTING OUTPUTS
    if (plot.display) {
      suppressWarnings(print(p3))
    }
    if (verbose >= 2) {
      # if (n.invariant > 0) {
        print(df)
      # } else {
      #   print(df[, c(
      #     "pop",
      #     "n.Ind",
      #     "n.Loc",
      #     "polyLoc",
      #     "monoLoc",
      #     "all_NALoc",
      #     "Ho",
      #     "HoSD",
      #     "He",
      #     "HeSD",
      #     "uHe",
      #     "uHeSD",
      #     "FIS",
      #     "FISSD"
      #   )], row.names = FALSE)
      # }
    }
  }
  
  ########### FOR METHOD BASED ON INDIVIDUAL
  
  if (method == "ind") {
    if (verbose >= 2) {
      cat(report("  Calculating observed heterozygosity for individuals\n"))
      cat(report(
        "  Note: No adjustment for invariant loci (n.invariant set to 0)\n"
      ))
    }
    # Convert to matrix
    m <- as.matrix(x)
    
    # For each individual: gametic heterozygosity d(k - d) / choose(k, 2)
    # averaged over its scored loci, and the shares of the two homozygotes
    # (dosage 0 and dosage k). k recycles down the columns, one per row.
    k <- as.numeric(ploidy(x))
    c.nloc <- unname(rowSums(!is.na(m)))
    c.hets <- unname(rowMeans(m * (k - m) / choose(k, 2), na.rm = TRUE))
    c.hom0 <- unname(rowMeans(m == 0, na.rm = TRUE))
    c.hom2 <- unname(rowMeans(m == k, na.rm = TRUE))

    # Join the sample sizes with the heterozygosities
    df <-
      cbind.data.frame(x@ind.names, c.hets, c.hom0, c.hom2,c.nloc)
    names(df) <-
      c("ind.name", "Ho", "f.hom.ref", "f.hom.alt","n.Loc")
    
    # Boxplot
    if (plot.display) {
      upper <- ceiling(max(df$Ho) * 10) / 10
      p1 <-
        ggplot(df, aes(y = Ho)) +
        geom_boxplot(color = plot.colors.ind[1], fill = plot.colors.ind[2]) +
        coord_flip() +
        plot.theme +
        xlim(range = c(-1, 1)) +
        ylim(0, upper) +
        ylab(" ") +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        ggtitle("Observed Heterozygosity by Individual")
      
      # Histogram
      p2 <-
        ggplot(df, aes(x = Ho)) +
        geom_histogram(bins = 25,
                       color = plot.colors.ind[1],
                       fill = plot.colors.ind[2]) +
        coord_cartesian(xlim = c(0, upper)) +
        xlab("Observed heterozygosity") +
        ylab("Count") +
        plot.theme
    }
    
    # Outliers are computed from the data (Tukey boxplot statistics, the
    # same rule ggplot uses) so the verbose >= 3 report does not depend
    # on the plot being displayed
    outliers_temp <- grDevices::boxplot.stats(df$Ho)$out
    outliers <-
      data.frame(ID = as.character(df$ind.name[df$Ho %in% outliers_temp]),
                 Ho = df$Ho[df$Ho %in% outliers_temp])

    # OUTPUT REPORT
    if (verbose >= 3) {
      cat("Reporting Heterozygosity by Individual\n")
      cat("No. of loci =", nLoc(x), "\n")
      cat("No. of individuals =", nInd(x), "\n")
      cat("  Minimum Observed Heterozygosity: ",
          round(min(df$Ho), 6),
          "\n")
      cat("  Maximum Observed Heterozygosity: ",
          round(max(df$Ho), 6),
          "\n")
      cat("  Average Observed Heterozygosity: ",
          round(mean(df$Ho), 6),
          "\n\n")
      cat("  Results returned as a dataframe\n\n")
      if (nrow(outliers) == 0) {
        cat("  No outliers detected\n\n")
      } else {
        cat("  Outliers detected\n")
        print(outliers)
        cat("\n")
      }
    }
    
    # PRINTING OUTPUTS
    if (plot.display) {
      p3 <- (p1 / p2) + plot_layout(heights = c(1, 4))
      print(p3)
    }
    if (verbose >= 2) {
      if(subsample.pop==TRUE){
        print(res_sub, row.names = FALSE)
      }
      print(df, row.names = FALSE)
    }
  }
  
  # Optionally save the plot ---------------------

  if(!is.null(plot.file)){
    if (exists("p3", inherits = FALSE)) {
      tmp <- utils.plot.save(p3,
                             dir=plot.dir,
                             file=plot.file,
                             verbose=verbose)
    } else if (verbose >= 1) {
      cat(warn(
        "  Warning: plot.file specified but no plot was built (plot.display is FALSE); nothing saved\n"
      ))
    }
  }
  
  if (verbose >= 3) {
    cat(report("  Returning a dataframe with heterozygosity values\n"))
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  if(subsample.pop==TRUE){
   return(invisible(list(subsample = res_sub, results = df)))
  }else{
  return(invisible(df))
  }

}
