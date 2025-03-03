#' Report SNP diversity from a genlight object, with reference to Ma, Z., Li, L., & Zhang, Y. P. (2020). Defining individual-level 
#' genetic diversity and similarity profiles. Scientific reports, 10(1), 5805.
#'
#' This function needs package adegenet, please install it. 
#' @param x A genlight file (works only for diploid data) [required].
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]..
#' @param plot.file Name for the RDS binary file to save (base name only, 
#' exclude extension) [default NULL].
#' @param level The types of SNP diversity to report. [default 'alpha', also accept 'beta', 'gamma'].
#' @param order The number of order to report. Starts from 0. [default 5].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' details
#' \itemize{
#' \item SNP diversity per individual}
#' @export
#' @author Ching Ching Lau (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @references 
#' \itemize{
#' \item Ma, Z., Li, L., & Zhang, Y. P. (2020). Defining individual-level 
#' genetic diversity and similarity profiles. Scientific reports, 10(1), 5805.}
#' @return A dataframe containing SNP diversity per individual
#' @examples
#' \dontrun{
#' obj <- gl.report.shannon(gl)
#' }

gl.report.shannon <- function(x,
                                plot.display = TRUE,
                                plot.theme = theme_dartR(),
                                plot.dir = NULL,
                                plot.file = NULL,
                                level = "alpha",
                                order=5,
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
  #if (is.null(pop(x)) |
  #    is.na(length(pop(x))) | length(pop(x)) <= 0) {
  #  if (verbose >= 2) {
  #    cat(
  #      warn(
  #        "  No population assignments detected,
  #                           individual assigned to a single population
  #                      labelled 'pop1'\n"
  #      )
  #    )
  #  }
  #  pop(x) <- array("pop1", dim = nInd(x))
  #  pop(x) <- as.factor(pop(x))
  #}
  
  # Codes of calculation are from Online Supplementary Information (OSI) for:
  # Ma ZS, Li LW and Zhang YP (2019) Defining Individual-Level Genetic
  # Diversity and Similarity Profiles. Scientific Reports
  
  d.chao <- function(A, lev, q) {
    tot <- sum(A)
    eA <- A / tot
    eA <- eA[eA > 0]
    if (is.vector(A)) {
      cA <- A
      N <- 1
    } else{
      cA <- colSums(A)
      N <- nrow(A)
    }
    ecA <- cA / tot
    ecA <- ecA[ecA > 0]
    if (lev == 'alpha') {
      if (q != 1) {
        Da <- (1 / N) * (sum(eA ^ q)) ^ (1 / (1 - q))
        D.value <- Da
      } else{
        Da <- exp(-sum(eA * log(eA)) - log(N))
        D.value <- Da
      }
    }
    if (lev == 'beta') {
      D.value <- d.chao(A, lev = 'gamma', q) / d.chao(A, lev = 'alpha', q)
    }
    if (lev == 'gamma') {
      if (q != 1) {
        Dg <- (sum(ecA ^ q)) ^ (1 / (1 - q))
        D.value <- Dg
      } else{
        Dg <- exp(-sum(ecA * log(ecA)))
        D.value <- Dg
      }
    }
    D.value
  }
  
  # SNP diversity
  if (verbose >= 2) {
    cat(
      report(
        "  Calculating SNP diversity for each individual\n"
      )
    )
  }
  
  x_mat <- as.matrix(x)
  ID <- rownames(x_mat)
  div_mat <- matrix(0, nrow(x_mat), order)
  list_order <- order-1
  for (n in 1:nrow(x_mat)) {
    otu <- x_mat[n, ]
    otu <- otu[otu > 0]
    otu <- otu[!is.na(otu)]
    for (q in 0:list_order) {
      div_mat[n, q + 1] <- d.chao(A = otu, lev = level, q)
    }
  }
  
  # output
  div_mat <- as.data.frame(div_mat)
  div_mat <- as.data.frame(cbind(ID, div_mat))
  colnames(div_mat) <- c("ID", paste0("q", c(0:list_order)))
  
  div_mat2 <- reshape2::melt(div_mat)
  Ord <- SNP_diversity <- NA
  colnames(div_mat2)[c(2,3)] <- c("Ord", "SNP_diversity")
  
  
  p1 <-
    ggplot(div_mat2, aes(
      x = ID,
    )) + geom_bar(aes(y = SNP_diversity, fill=Ord), position = "dodge",
                  stat = "identity") + plot.theme + theme(
                    axis.ticks.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.ticks.y = element_blank(),
                    legend.position = "none"
                  ) + facet_grid(~Ord) +
    labs(fill = "Order", y=paste0(level, " diversity")) +
    ggtitle(paste0(level, " diversity per individual for different orders"))

  
  # Optionally save the plot ---------------------
  
  if(!is.null(plot.file)){
    tmp <- utils.plot.save(p1,
                           dir=plot.dir,
                           file=plot.file,
                           verbose=verbose)
  }
  
  if (verbose >= 3) {
    cat(report("  Returning a matrix with SNP diversity values\n"))}
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # PRINTING OUTPUTS
  if (plot.display) {
    suppressWarnings(print(p1))}
  
  # RETURN
  return(invisible(div_mat))}
