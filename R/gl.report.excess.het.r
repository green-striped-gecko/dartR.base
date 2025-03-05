#' @name gl.report.excess.het
#' @title Report loci with excess of heterozygosity
#' @family unmatched report
#'
#' @description Calculates excess of heterozygosity in a genlight object

#' @param x Name of the genlight object containing the SNP data [required].
#' @param Yates Boolean for Yates's continuity correction. [default FALSE]
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors Vector with two color names for the borders and fill
#' [default c("#2171B5", "#6BAED6")].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()].
#' @param plot.file Name for the RDS binary file to save (base name only, 
#' exclude extension) [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#'
#' @details
#' Loci with observed heterozygosity larger than 0.5 and expected heterozygosity would be indicated as excess (p-value <= 0.05). 
#' You can remove the loci with excess of heterozygosity from genlight object using \code{\link{gl.filter.excess.het}}
#'
#'\strong{ Function's output }
#' 
#' If a plot.file is given, the ggplot arising from this function is saved as an 
#' "RDS" binary file using saveRDS(); can be reloaded with readRDS(). A file 
#' name must be specified for the plot to be saved.
#' 
#' If a plot directory (plot.dir) is specified, the ggplot binary is saved to 
#' that directory; otherwise to the tempdir(). 
#'  
#'  Examples of other themes that can be used can be consulted in: 
#'  \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'  A color vector can be obtained with gl.select.colors() and then passed to
#'  the function with the plot.colors parameter.
#' @author Jesús Castrejón-Figueroa, Diana A Robledo-Ruiz (Custodian: Ching Ching Lau) -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @references
#' \itemize{
#' \item https://github.com/drobledoruiz/conservation_genomics/tree/main/filter.excess.het
#' \item Robledo‐Ruiz, D. A., Austin, L., Amos, J. N., Castrejón‐Figueroa, J., Harley, D. K., Magrath, M. J., Sunnucks, P. & Pavlova, A. (2023). 
#' Easy‐to‐use R functions to separate reduced‐representation genomic datasets into sex‐linked and autosomal loci, 
#' and conduct sex assignment. Molecular Ecology Resources.
#' }
#'
#' @examples
#' filtered.table <- gl.report.excess.het(x = LBP, Yates = TRUE)
#' @seealso \code{\link{gl.filter.excess.het}}
#' @importFrom stats aggregate
#' @export
#' @return 1. Table with information of excessively-heterozygous loci \cr
#' 2. Two plots of heterozygosity of the loci before and after filtering (i.e. without excessively heterozygous loci).\cr
#' 3. A vector with the names of loci to be remove by \code{\link{gl.filter.excess.het}}

gl.report.excess.het <- function(x,
                                 Yates=FALSE,
                                 plot.display = TRUE,
                                 plot.theme = theme_dartR(),
                                 plot.colors = NULL,
                                 plot.file = NULL,
                                 plot.dir = NULL,
                                 verbose = NULL) {
  
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  if(verbose==0){plot.display <- FALSE}
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir,verbose=0)
  
  # SET COLOURS
  if (is.null(plot.colors)) {
    plot.colors <- c("#2171B5", "#6BAED6")
  } else {
    if (length(plot.colors) > 2) {
      if (verbose >= 2) {
        cat(warn("  More than 2 colors specified, only the first 2 are used\n"))
      }
      plot.colors <- plot.colors[1:2]
    }
  }
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  
  # DO THE JOB
  
  if(!Yates) {
    cc = 0
  } else {
    cc = 0.5
  }
  
  # Start BEFORE plot and results table with observed data
  # Plot
  if (verbose >= 2) {
    cat(report("  Building BEFORE-filtering plot"))
  }
  
  gen <- as.data.frame(t(as.matrix(x)))
  n0 <- rowSums(gen == 0, na.rm = TRUE)
  n1 <- rowSums(gen == 1, na.rm = TRUE)
  n2 <- rowSums(gen == 2, na.rm = TRUE)
  
  fhe <- n1/(n0 + n1 + n2)
  Index <- loc <-NULL
  p1 <- ggplot(data.frame(loc=fhe,Index=1:nLoc(x)), aes(x= Index,y = loc)) + 
  geom_point(fill=plot.colors[2], colour=plot.colors[1], pch=21) +
  plot.theme +
  ylim(0,1) + 
  theme(axis.ticks.x = element_blank()) +
    labs(title="BEFORE", x="Index", y="Locus heterozygosity")

  # Results per population
  populations <- as.vector(unique(x@other$ind.metrics$pop))
  n0   <- vector()
  n1   <- vector()
  n2   <- vector()
  pops <- vector()
  loci <- vector()
  
  for (pop in populations) {
    pop.gl <- x[x@other$ind.metrics$pop == pop,]
    gen    <- as.data.frame(t(as.matrix(pop.gl)))
    n0     <- c(n0, rowSums(gen == 0, na.rm = TRUE))
    n1     <- c(n1, rowSums(gen == 1, na.rm = TRUE))
    n2     <- c(n2, rowSums(gen == 2, na.rm = TRUE))
    loci   <- c(loci, rownames(gen))
    pops   <- c(pops, rep(pop, dim(gen)[[1]]))
  }
  
  table <- data.frame(
    loci = loci,
    pop  = pops,
    n0   = n0,
    n1   = n1,
    n2   = n2,
    Hobs = n1 / (n0 + n1 + n2)
  )
  
  # Function to test HW with chisq
  HWchsq <- function(n0, n1, n2, cc) {
    n    <- n0 + n1 + n2
    p    <- (2 * n0 + n1) / (2 * n)
    q    <- 1 - p
    en0  <- n * p * p
    en1  <- 2 * n * p * q
    en2  <- n * q * q
    Hexp <- en1 / (en0 + en1 + en2)
    chsq <-
      (abs(n0 - en0) - cc) ** 2 / en0 + (abs(n1 - en1) - cc) ** 2 / en1 + (abs(n2 - en2) - cc) **
      2 / en2
    return(c(en0, en1, en2, Hexp, chsq))
  }
  
  # Function to test HW for each locus
  HW.chsqtest <- function(table, cc = cc) {
    table$En0     <- NA
    table$En1     <- NA
    table$En2     <- NA
    table$Hexp    <- NA
    table$chsq    <- NA
    table$p.value <- NA
    
    for (i in 1:nrow(table)) {
      n0 <- table$n0[i]
      n1 <- table$n1[i]
      n2 <- table$n2[i]
      
      chsq <- HWchsq(n0, n1, n2, cc)
      
      table[i, 'En0']     <- chsq[1]
      table[i, 'En1']     <- chsq[2]
      table[i, 'En2']     <- chsq[3]
      table[i, 'Hexp']    <- chsq[4]
      table[i, 'chsq']    <- chsq[5]
      table[i, 'p.value'] <- pchisq(chsq[5], 1, lower.tail = FALSE)
    }
    return(table)
  }
  
  # Apply chsq test only to loci with Het > 50%
  if (verbose >= 2) {
    cat(report("  Testing loci for high heterozygosity..."))
  }
  
  table.filter <- na.omit(table[table$Hobs >= 0.5,])
  table.filter <- HW.chsqtest(table.filter, cc = cc)
  
  # Adjust p-values
  table.filter$p.adjusted <-
    p.adjust(table.filter$p.value, method = "fdr")
  
  # Keep in table only loci p.adj <= 0.05 and ObsHet >= ExpHet
  table.filter <-
    table.filter[table.filter$p.adjusted <= 0.05 &
                   table.filter$n1 >= table.filter$En1, ]
  
  # Check existence of excessively heterozygous loci
  if (nrow(table.filter) == 0){  
    message("No excessively-heterozygous loci found.")  
  } else {
    rownames(table.filter) <- 1:nrow(table.filter)}
    
    # Remove highly-het loci from new filtered gl
    gl.filter <- x[, !(x$loc.names %in% table.filter$loci )]
    
    # Remove highly-het loci from new filtered gl
    x2 <- x[,!(x$loc.names %in% table.filter$loci)]
    
    # Remove highly-het loci from new gl loci metadata
    x2@other$loc.metrics <-
      x@other$loc.metrics[!(x$loc.names %in% table.filter$loci),]
    
    # AFTER plot with filtered gl
    if (verbose >= 2) {
      cat(report("  Building AFTER-filtering plot"))
    }
    
    gen <- as.data.frame(t(as.matrix(x2)))
    n0 <- rowSums(gen == 0, na.rm = TRUE)
    n1 <- rowSums(gen == 1, na.rm = TRUE)
    n2 <- rowSums(gen == 2, na.rm = TRUE)
    
    fhe <- n1 / (n0 + n1 + n2)

      
    p2 <- ggplot(data.frame(loc=fhe,Index=1:nLoc(x2)), aes(x= Index,y = loc)) + 
      geom_point(fill=plot.colors[2], colour=plot.colors[1], pch=21) +
      plot.theme +
      ylim(0,1) + 
      theme(axis.ticks.x = element_blank()) +
      labs(title="AFTER", x="Index", y="Locus heterozygosity")
    
    # using package patchwork
    p3 <- (p1 / p2) 
    if (plot.display) {print(p3)}
    
    if (verbose >= 3) {
        print(list('results.table' = table.filter,
                   'removed.loci'  = unique(table.filter$loci)))
      }
  
    ################## 6. Output

    cat(report("\n  **FINISHED** Detected ", x@n.loc - gl.filter@n.loc, " excessively-heterozygous loci."))
    
  # Optionally save the plot ---------------------
  
  if(!is.null(plot.file)){
    tmp <- utils.plot.save(p3,
                           dir=plot.dir,
                           file=plot.file,
                           verbose=verbose)
  }
  
  if (verbose >= 3) {
    cat(report("  Returning a dataframe with loci heterozygosity values\n"))
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("\nCompleted:", funname, "\n"))
  }
  
  # RETURN
  return(invisible(list('results.table' = table.filter,
                        'removed.loci'  = unique(table.filter$loci))))
} 
