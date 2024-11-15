#' Report allelic richness per population from a genlight object
#'
#' This function needs package adegenet, please install it. 
#' @param x A genlight file (works only for diploid data) [required].
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()].
#' @param plot.file Name for the RDS binary file to save (base name only, 
#' exclude extension) [default NULL].
#' @param error.bar Statistic to be plotted as error bar either "SD" (standard 
#' deviation) or "SE" (standard error) or "CI" (confidence intervals)
#'  [default "SD"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' details
#' \itemize{
#' \item raw allele count
#' \item allelic richness}
#' @export
#' @author Ching Ching Lau (Post to \url{https://groups.google.com/d/forum/dartr})
#' @references 
#' \itemize{
#' \item El Mousadik, A., & Petit, R. J. (1996). High level of genetic differentiation for allelic richness among populations of the argan tree 
#' [Argania spinosa (L.) Skeels] endemic to Morocco. Theoretical and applied genetics, 92, 832-839.}
#' @return A dataframe containing richness per site, richness per population, raw reference allele count,
#' raw alternate allele count.
#' @examples
#'  gl.report.allelerich(possums.gl)

gl.report.allelerich <- function(x,
                                plot.display = TRUE,
                                plot.theme = theme_dartR(),
                                plot.dir = NULL,
                                plot.file = NULL,
                                error.bar = "SD",
                                verbose = 2) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  pkg <- c("dplyr", "tidyr")
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  } 
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir,verbose=0)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jackson",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <-
    utils.check.datatype(x, accept = "SNP", verbose = verbose)
  
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
  sgl <- seppop(x)
  
  min_pop <- NULL
  allele_count_all <- NULL
  
  for (l in 1:length(sgl)) {
    # convert genlight to SNP matrix
    m <- as.matrix(sgl[[l]])
    allele_count <- reshape2::melt(m) 
    allele_count$pop <- sgl[[l]]$pop[1]
    colnames(allele_count)[c(2,3)] <- c("site", "genotype")
    
    # summarise allele count and richness
    allele_count2 <-allele_count[,-1] %>% group_by(site, genotype, pop) %>% count()
    allele_count2 <- na.omit(allele_count2)
    colnames(allele_count2)[4] <- "n"
    allele_count2$ref_allele <- NA
    allele_count2$alt_allele <- NA
    # calculate number of alternate and reference allele, assuming diploid
    allele_count2[which(allele_count2$genotype==0),'ref_allele'] <- allele_count2[which(allele_count2$genotype==0),'n']*2
    allele_count2[which(allele_count2$genotype==1),'ref_allele'] <- allele_count2[which(allele_count2$genotype==1),'n']
    allele_count2[which(allele_count2$genotype==2),'ref_allele'] <- 0
    allele_count2[which(allele_count2$genotype==0),'alt_allele'] <- 0
    allele_count2[which(allele_count2$genotype==1),'alt_allele'] <- allele_count2[which(allele_count2$genotype==1),'n']
    allele_count2[which(allele_count2$genotype==2),'alt_allele'] <- (allele_count2[which(allele_count2$genotype==2),'n']*2)
    allele_count3 <- allele_count2 %>% dplyr::group_by(site) %>% summarise(pop=pop, all_ref_allele=sum(ref_allele), all_alt_allele=sum(alt_allele))
    allele_count_all <- rbind(allele_count_all, distinct(allele_count3, site, .keep_all = T) )
    popsize_per_snp <- allele_count_all %>% group_by(site, pop) %>% summarise(raw_count=(all_ref_allele + all_alt_allele))
    min_pop <-  min(popsize_per_snp %>% group_by(pop) %>% summarise(min_pop=min(raw_count))%>%select(min_pop)) }
  
  
  #min_pop <- min(popsize_per_snp$raw_count)
  
  # calculate allele richness
  summary_pop_allele <- ((merge(popsize_per_snp, allele_count_all)))
  summary_pop_allele$ref_site_richness <- (1-choose(((summary_pop_allele$raw_count*2)-(summary_pop_allele$all_alt_allele*2)),min_pop)/choose((summary_pop_allele$raw_count*2),min_pop))
  summary_pop_allele$alt_site_richness <- (1-choose(((summary_pop_allele$raw_count*2)-(summary_pop_allele$all_ref_allele*2)),min_pop)/choose((summary_pop_allele$raw_count*2),min_pop))
  summary_pop_allele$`sum_site_richness` <- (summary_pop_allele$ref_site_richness+summary_pop_allele$alt_site_richness)
  #summary_pop_allele=NULL

  # output
  richness <- tidyr::pivot_wider(summary_pop_allele[,c('pop','sum_site_richness', 'site')], names_from=pop, values_from=sum_site_richness)
  richness_summary <- data.frame(sum_richness=colSums(richness[,-c(1)], na.rm = T), mean_richness=colMeans(richness[,-c(1)], na.rm = T))
  richness_summary$pop <- rownames(richness_summary)
  
  # assign original populaiton size
  richness_summary$popsize <- NA
  for (i in 1:nrow(richness_summary)){
    richness_summary$popsize[i] <- nInd(sgl[[richness_summary$pop[i]]])
  }
  raw_count_ref <- tidyr::pivot_wider(summary_pop_allele[,c('pop','all_ref_allele', 'site')], names_from=pop, values_from=all_ref_allele)
  raw_count_alt <- tidyr::pivot_wider(summary_pop_allele[,c('pop','all_alt_allele', 'site')], names_from=pop, values_from=all_alt_allele)
  res <- list(richness, richness_summary, raw_count_ref, raw_count_alt)
  names(res) <- c("Richness per site", "Richness per population", "Raw reference allele count", "Raw alternate allele count")
  
    # error bar
    if(error.bar=="SD") {
      pop_list_plot_error <- apply(richness[,-c(1)],2, sd, na.rm=T)
    }
  
    if(error.bar=="SE") {
      len_of_m <- apply(richness[,-c(1)], 2, length)
      pop_list_plot_error <- apply(richness[,-c(1)],2, sd, na.rm=T)/sqrt(len_of_m)
    }
    
  
    p1 <-
      ggplot(richness_summary, aes(
        x = pop,
        y = sum_richness,
        fill = pop
      )) + geom_bar(position = "dodge",
                  stat = "identity",
                  color = "black") + plot.theme + theme(
                    axis.ticks.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position = "none"
                  ) +
      labs(fill = "Population") +
     ggtitle("Sum allelic richness by Population")
  
    p2 <-
      ggplot(richness_summary, aes(
        x = paste(pop,"n=", popsize),
        y = mean_richness,
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
      labs(fill = "Population") +
      ggtitle("Mean allelic richness by Population")
  
    if(error.bar=="SD"){
      p2 <- p2 + 
        geom_errorbar(aes(ymin = mean_richness - pop_list_plot_error, 
                        ymax = mean_richness + pop_list_plot_error), 
                    width=0.5)+
        ggtitle(label = "Mean allelic richness by Population",
              subtitle = "Error bars show Standard Deviation")
    }
  
    if(error.bar=="SE"){
      p2 <- p2 + 
        geom_errorbar(aes(ymin = mean_richness - pop_list_plot_error, 
                        ymax = mean_richness + pop_list_plot_error), 
                    width=0.5) +
        ggtitle(label = "Mean allelic richness by Population",
              subtitle = "Error bars show Standard Error")
    }
  
    p3 <- (p1 / p2)
    
# Optionally save the plot ---------------------
    
    if(!is.null(plot.file)){
      tmp <- utils.plot.save(p3,
                             dir=plot.dir,
                             file=plot.file,
                             verbose=verbose)
    }
    
    if (verbose >= 3) {
      cat(report("  Returning a dataframe with allelic richness values\n"))}
      
    # FLAG SCRIPT END
    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }
  
  # PRINTING OUTPUTS
    if (plot.display) {
      suppressWarnings(print(p3))}
    
  # RETURN
      return(invisible(res))}
