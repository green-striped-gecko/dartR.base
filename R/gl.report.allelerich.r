#' @name gl.report.allelerich
#' @title Reports allelic richness per population from a genlight object
#'
#' @description
#' This function calculates allelic richness across populations for SNP data 
#' in a genlight object, using a rarefaction-based approach 
#' adapted from El Mousadik and Petit (1996). By standardizing the expected number 
#' of distinct alleles to a fixed sample size, it allows meaningful comparisons 
#' among populations with different sampling depths. Because allelic richness 
#' captures the contribution of rare alleles more effectively than heterozygosity-based 
#' indices, it can better reveal subtle patterns of genetic diversity, especially 
#' in conservation contexts where low-frequency alleles may hold particular importance. 
#' Optional bootstrapping is provided to generate measures of uncertainty such as 
#' confidence intervals, standard error, or standard deviation.
#'
#' @family unmatched report
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param nboots Number of bootstrap replicates to obtain confidence intervals
#'   [default 0].
#' @param boot.method Character specifying the bootstrap strategy: "ind" to
#'   resample individuals or "loc" to resample loci [default "loc"].
#' @param conf Numeric specifying the confidence level for the interval 
#'   estimation [default 0.95].
#' @param CI.type Character specifying the type of nonparametric confidence 
#'   interval: "norm", "basic", "perc", or "bca" [default "bca"].
#' @param plot.display Logical indicating if a plot should be produced 
#'   [default TRUE].
#' @param plot.theme A \pkg{ggplot2} theme to style the plots 
#'   [default \code{theme_dartR()}].
#' @param plot.colors.pop A color palette for population plots or a list with
#' as many colors as there are populations in the dataset 
#' [default gl.colors("dis")].
#' @param plot.dir Directory in which to save plot RDS files [default is the 
#'   current working directory or \code{tempdir()}].
#' @param plot.file Base name (excluding extension) for the RDS file to save 
#'   [default NULL].
#' @param error.bar Character specifying the statistic used for error bars: 
#'   "SD" (standard deviation), "SE" (standard error), or "CI" (confidence 
#'   intervals) [default "SD"].
#' @param verbose Numeric controlling the level of printed messages: 
#'   0 = silent or fatal errors; 1 = begin/end; 2 = progress log; 3 = progress 
#'   and summary; 5 = full report. Defaults to 2 unless changed via 
#'   \code{gl.set.verbosity}.
#'
#' @details
#' \strong{Allelic Richness via Rarefaction:}
#' 
#' This function applies a rarefaction technique to standardize allelic 
#' richness. This approach mitigates biases that arise from differences 
#' in sampling depth and highlights the presence of less frequent alleles, 
#' which can be disproportionately important for long-term evolutionary 
#' potential and conservation. Rarefaction essentially calculates the 
#' expected number of distinct alleles you would observe in a subsample 
#' of fixed size.
#' 
#' \strong{Interpretation in Conservation Genetics:}
#' 
#' Because allelic richness emphasizes rare alleles, it often reveals 
#' patterns of genetic differentiation that may be underestimated by 
#' measures relying primarily on allele-frequency distributions. In studies 
#' aiming to preserve or identify unique variants, allelic richness can 
#' highlight the most genetically distinct populations or pinpoint those 
#' at greater risk of losing rare allelic variants.
#' 
#'   \strong{Error bars}
#'  
#'  The best method for presenting or assessing genetic statistics depends on 
#'  the type of data you have and the specific questions you're trying to 
#'  answer. Here's a brief overview of when you might use each method:
#'  
#'   \strong{1. Confidence Intervals ("CI"):}
#'   
#'- Usage: Often used to convey the precision of an estimate.
#'  
#'- Advantage: Confidence intervals give a range in which the true parameter 
#'  (like a population mean) is likely to fall, given the data and a specified 
#'  probability (like 95\%).
#'  
#'- In Context: For genetic statistics, if you're estimating a parameter,
#' a 95\% CI gives you a range in which you're 95\% confident the true parameter
#'    lies.
#'  
#'   \strong{2. Standard Deviation ("SD"):}
#'   
#'- Usage: Describes the amount of variation from the average in a set of data.
#'  
#'- Advantage: Allows for an understanding of the spread of individual data
#'   points around the mean.
#'   
#'- In Context: If you're looking at the distribution of a quantitative trait 
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
#'   In practice, geneticists often use a combination of these methods to 
#'   analyze and present their data, depending on their research questions and 
#'   the nature of the data.
#'   
#'  \strong{Confidence Intervals}
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
#' The studentised bootstrap interval ("stud") was not included in the CI types
#'  because it is computationally intensive, it may produce estimates outside
#'  the range of plausible values and it has been found to be erratic in
#'  practice, see for example the "Studentised (t) Intervals" section in:
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
#' \strong{Confidence intervals that do not encompass the mean}
#'  
#' When using bootstrap methods to estimate confidence intervals for rarefied
#'  metrics, the resampling distribution can be skewed—especially if the sample
#'   size is small. Skewness means the bootstrap‐estimated distribution may 
#'   shift away from the sample mean, so percentile (or BCa) confidence 
#'   intervals may legitimately exclude that mean. This does not necessarily 
#'   indicate an error but can reflect genuine asymmetry in the data.
#'
#' @return
#' A list of data frames containing:
#' \itemize{
#'   \item Allelic Richness per site (corrected by rarefaction).
#'   \item Allelic Richness per population (summarized across sites).
#'   \item Raw reference allele counts.
#'   \item Raw alternate allele counts.
#' }
#'
#' @references
#' El Mousadik, A. & Petit, R. J. (1996). High level of genetic 
#' differentiation for allelic richness among populations of the argan tree 
#' [\emph{Argania spinosa} (L.) Skeels] endemic to Morocco. 
#' \emph{Theoretical and Applied Genetics}, 92, 832–839.
#'
#' @author Author(s): Ching Ching Lau -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @export
#'
#' @importFrom dplyr group_by select summarise distinct
#'
#' @examples
#'   # Example usage:
#'   if (isTRUE(getOption("dartR_fbm"))) possums.gl <- gl.gen2fbm(possums.gl)
#'   results <- gl.report.allelerich(possums.gl)
#'   print(results)
#'

gl.report.allelerich <- function(x,
                                 nboots = 0,
                                 boot.method = "loc",
                                 conf = 0.95,
                                 CI.type = "bca",
                                 plot.display = TRUE,
                                 plot.theme = theme_dartR(),
                                 plot.colors.pop = gl.colors("dis"),
                                 plot.dir = NULL,
                                 plot.file = NULL,
                                 error.bar = "SD",
                                 verbose = NULL) {
  
  # # setting parallel
  # if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
  #   parallel <- "multicore"
  # }
  # ## if windows
  # if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
  #   parallel <- "snow"
  # }
  min_sample <- overall_min <- r_ref <- r_alt <- r_total <- 
    corrected_richness <- sum_corrected_richness <- 
    mean_corrected_richness <- NULL
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  pkg <- c("dplyr", "tidyr","reshape2")
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jackson",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)
  
  # check population assignment
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
  
  if( nboots == 0 & error.bar == "CI" ){
    cat(error(
      "  Number of boostraps ('nboots' parameter) must be > 0 to calculate confidence 
   intervals \n"))
    stop()
  }
  
  if( nboots > 0){
  error.bar <- "CI"
  }
  
  #define some global variables...
  site <- genotype <- ref_allele <- alt_allele <- raw_count <-  
    sum_site_richness <- sum_richness <- popsize <- mean_richness <- 
    all_ref_allele <- all_alt_allele <-  NA
  
  # ALLELIC RICHNESS
  if (verbose >= 2) {
    cat(
      report(
        "  Calculating Allelic Richness, averaged across
                    loci, for each population\n"
      )
    )
  }
  
  # Split the genlight object into a list of populations
  pop_list <- seppop(x)

  allele_site_summary <- lapply(names(pop_list), function(pop_name) {
    pop_data <- pop_list[[pop_name]]
    # Convert genlight object to a matrix and reshape to long format
    m <- as.matrix(pop_data)
    allele_df <- reshape2::melt(m, varnames = c("ind", "site"), value.name = "genotype") %>%
      dplyr::mutate(pop = pop_name)
    
    # Summarise counts per genotype for each SNP and compute allele counts:
    # - Genotype 0: 2 copies of the reference allele, 0 copies of the alternate allele.
    # - Genotype 1: 1 copy each of reference and alternate.
    # - Genotype 2: 0 copies of the reference allele, 2 copies of the alternate allele.
    allele_summary <- allele_df %>%
      group_by(site, genotype, pop) %>%
      dplyr::tally(name = "n") %>%
      stats::na.omit() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        ref_allele = dplyr::case_when(
          genotype == 0 ~ n * 2,
          genotype == 1 ~ n,
          genotype == 2 ~ 0
        ),
        alt_allele = dplyr::case_when(
          genotype == 0 ~ 0,
          genotype == 1 ~ n,
          genotype == 2 ~ n * 2
        )
      )
    
    # Summarise total allele counts per site for the population
    site_summary <- allele_summary %>%
      dplyr::group_by(site, pop) %>%
      summarise(
        all_ref_allele = sum(ref_allele,na.rm = TRUE),
        all_alt_allele = sum(alt_allele,na.rm = TRUE),
        raw_count = sum(ref_allele, alt_allele,na.rm = TRUE),
        .groups = "drop"
      )
    
    return(site_summary)
  })
  
  # Combine data from all populations into one data frame
  allele_count_all <- dplyr::bind_rows(allele_site_summary)
  
  # Determine the overall minimum sample size (i.e., the smallest total number
  # of allele copies) across all sites and populations. This minimum is used as 
  # the subsample size (n) in the rarefaction formula.
  min_pop <- allele_count_all %>%
    dplyr::group_by(pop) %>%
    summarise(min_sample = min(raw_count), .groups = "drop") %>%
    summarise(overall_min = min(min_sample)) %>%
    dplyr::pull(overall_min)
  
  # ---------------------------------------------------------------------------
  # Integration of Rarefaction Equations:
  #
  # The rarefaction method allows for the standardization of allelic richness
  # across samples of different sizes. Here we incorporate the equations from the paper:
  #
  # For each site:
  #   - Let N be the total number of gene copies (raw_count),
  #   - N_ref be the count of the reference allele (all_ref_allele),
  #   - N_alt be the count of the alternate allele (all_alt_allele),
  #   - n (here min_pop) be the subsample size.
  #
  # Then, for each allele type, the expected probability of occurrence in a subsample
  # is computed as:
  #
  #   r_ref = 1 - choose(N - N_ref, n) / choose(N, n)
  #   r_alt = 1 - choose(N - N_alt, n) / choose(N, n)
  #
  # The expected total number of alleles in the subsample is:
  #
  #   r(n) = r_ref + r_alt
  #

  allele_richness <- allele_count_all %>%
    dplyr::mutate(
      r_ref = 1 - choose(raw_count - all_ref_allele, min_pop) / 
        choose(raw_count, min_pop),
      r_alt = 1 - choose(raw_count - all_alt_allele, min_pop) / 
        choose(raw_count, min_pop),
      r_total = r_ref + r_alt,
      corrected_richness = r_total
    )
  
  # Pivot the data to obtain corrected richness per site in wide format by
  # population.
  richness_per_site <- allele_richness %>%
    dplyr::select(pop, site, corrected_richness) %>%
    tidyr::pivot_wider(names_from = pop, values_from = corrected_richness)
  
  # Summarise the corrected richness per population (averaging over sites)
  richness_summary <- allele_richness %>%
    dplyr::group_by(pop) %>%
    summarise(
      sum_corrected_richness = sum(corrected_richness, na.rm = TRUE),
      mean_corrected_richness = mean(corrected_richness, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      tibble::tibble(pop = names(pop_list),
             popsize = sapply(pop_list, nInd)),
      by = "pop"
    )
  
  # Create wide-format tables for raw allele counts per site for reference and 
  # alternate alleles
  raw_count_ref <- allele_richness %>%
    dplyr::select(pop, site, all_ref_allele) %>%
    tidyr::pivot_wider(names_from = pop, values_from = all_ref_allele)
  
  raw_count_alt <- allele_richness %>%
    dplyr::select(pop, site, all_alt_allele) %>%
    tidyr::pivot_wider(names_from = pop, values_from = all_alt_allele)
  
  # The final output is a list containing:
  # - "Corrected Richness per site": Allelic richness per site computed by 
  # rarefaction.
  # - "Corrected Richness per population": Summaries of corrected richness 
  # averaged over sites.
  # - "Raw reference allele count": The original count of the reference allele.
  # - "Raw alternate allele count": The original count of the alternate allele.
  result <- list(
    "Allelic Richness per site" = richness_per_site,
    "Allelic Richness per population" = richness_summary,
    "Raw reference allele count" = raw_count_ref,
    "Raw alternate allele count" = raw_count_alt
  )
  
  # Round numeric columns in each data frame of the result list to 4 decimal
  # places
  result <- lapply(result, function(df) {
    df <- df %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., digits = 6)))
    as.data.frame(df)
  })
  
  npops <- nPop(x)
  # bootstrapping
  if (nboots > 0) {
    # Split the genlight object into a list of populations
    sgl <- seppop(x)
    pop_boot <- lapply(sgl, function(y) {
      df <- as.matrix(y)
      
      if(boot.method == "loc"){
        df <- t(df)
      }
      
      res_boots <- boot::boot(
        data = df,
        statistic = utils.allelic.richness,
        boot_method = boot.method,
        R = nboots,
        parallel = "no"
      )
      return(res_boots)
    })
    
    # confidence intervals
    
    # creating matrices to store CI
    nparams <- ncol(pop_boot[[1]]$t)
    res_CI <- replicate(npops,
                        as.data.frame(matrix(nrow = nparams, ncol = 2)),
                        simplify = FALSE)
    
    for (pop_n in 1:length(sgl)) {
      for (stat_n in seq_len(nparams)) {
        res_CI_tmp <- boot::boot.ci(
          boot.out = pop_boot[[pop_n]],
          conf = conf,
          type = CI.type,
          index = stat_n,
          t0 =  pop_boot[[pop_n]]$t0,
          t = pop_boot[[pop_n]]$t[, stat_n]
        )
        
        res_CI[[pop_n]][stat_n,] <-
          tail(as.vector(res_CI_tmp[[4]]), 2)
        
      }
    }
    
    res_CI <- Reduce(rbind,res_CI)
    colnames(res_CI) <- c("LCI", "HCI")
    
  }
  
  # error bar
  if (error.bar == "SD") {
    pop_list_plot_error <- apply(richness_per_site[, -c(1)], 2, sd, na.rm = T)
    result$`Allelic Richness per population`$SD <- pop_list_plot_error
    max_val <- max(result$`Allelic Richness per population`$SD + 
                     result$`Allelic Richness per population`$mean_corrected_richness)
  }
  
  if (error.bar == "SE") {
    len_of_m <- apply(richness_per_site[, -c(1)], 2, length)
    pop_list_plot_error <- apply(richness_per_site[, -c(1)], 2, sd, na.rm = T) /
      sqrt(len_of_m)
    result$`Allelic Richness per population`$SE <- pop_list_plot_error
    max_val <- max(result$`Allelic Richness per population`$SE + 
                     result$`Allelic Richness per population`$mean_corrected_richness)
  }
  
  if (error.bar == "CI") {
    result$`Allelic Richness per population` <- 
      cbind(result$`Allelic Richness per population`,
            res_CI)
    max_val <- max(res_CI$HCI)
  }
  
  # printing plots and reports assigning colors to populations
  if (is(plot.colors.pop, "function")) {
    colors_pops <- plot.colors.pop(length(levels(pop(x))))
  }
  
  if (!is(plot.colors.pop, "function")) {
    colors_pops <- plot.colors.pop
  }
    
  p1 <-
    ggplot(result[[2]], 
           aes(x = pop, y = sum_corrected_richness, fill = pop)) +
    geom_bar(position = "dodge",
             stat = "identity",
             color = "black") +
    plot.theme +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none" ) +
    labs(fill = "Population") +
    ggtitle("Sum allelic richness by Population")

  p2 <-
    ggplot(result[[2]], aes(
      x = paste(pop, "n=", popsize),
      y = mean_corrected_richness,
      fill = pop
    )) + geom_bar(position = "dodge",
                  stat = "identity",
                  color = "black") + plot.theme + theme(
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
    scale_fill_manual(values = colors_pops) +
    coord_cartesian(ylim = c(0.8,max_val)) + 
    labs(fill = "Population") +
    ggtitle("Mean allelic richness by Population")

  if (error.bar == "SD") {
    p2 <- p2 +
      geom_errorbar(
        aes(
          ymin = mean_corrected_richness - pop_list_plot_error,
          ymax = mean_corrected_richness + pop_list_plot_error
        ),
        width = 0.5
      ) +
      ggtitle(label = "Mean allelic richness by Population",
              subtitle = "Error bars show Standard Deviation")
  }

  if (error.bar == "SE") {
    p2 <- p2 +
      geom_errorbar(
        aes(
          ymin = mean_corrected_richness - pop_list_plot_error,
          ymax = mean_corrected_richness + pop_list_plot_error
        ),
        width = 0.5
      ) +
      ggtitle(label = "Mean allelic richness by Population",
              subtitle = "Error bars show Standard Error")
  }
  
  if (error.bar == "CI") {
    p2 <- p2 +
      geom_errorbar(
        aes(
          ymin = res_CI$LCI,
          ymax = res_CI$HCI
        ),
        width = 0.5
      ) +
      ggtitle(label = "Mean allelic richness by Population",
              subtitle = "Error bars show Confidence Intervals")
  }

  # p3 <- (p1 / p2)
  p3 <-  p2
  

  # Optionally save the plot ---------------------

  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p3,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }

  if (verbose >= 3) {
    cat(report("  Returning a dataframe with allelic richness values\n"))
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # PRINTING OUTPUTS
  if (plot.display) {
    suppressWarnings(print(p3))
  }
  
  # RETURN
  return(invisible(result))
}


