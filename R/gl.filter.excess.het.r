#' @name gl.filter.excess.het
#' @title Filters excessively-heterozygous loci from a genlight object
#' @family matched report
#'
#' @description Calculates excess of heterozygosity in a genlight object and remove those loci
#' @param x A genlight object containing the SNP genotypes [required].
#' @param Yates Whether to use Yates's continuity correction [default FALSE].
#' @param mono.rm Remove monomorphic loci after analysis is complete
#' [default FALSE].
#' @param recalc Recalculate the locus metadata statistics if any individuals
#' are deleted in the filtering [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
# Plot themes can be obtained from:
#  \itemize{
#  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#  }
#'
#'
#' @author Author(s): Jesús Castrejón-Figueroa, Diana A Robledo-Ruiz (Custodian: Ching Ching Lau) -- Post
#' to \url{https://groups.google.com/d/forum/dartr}
#' @references
#' \itemize{
#' \item https://github.com/drobledoruiz/conservation_genomics/tree/main/filter.excess.het
#' \item Robledo‐Ruiz, D. A., Austin, L., Amos, J. N., Castrejón‐Figueroa, J., Harley, D. K., Magrath, M. J., ... & Pavlova, A. (2023). 
#' Easy‐to‐use R functions to separate reduced‐representation genomic datasets into sex‐linked and autosomal loci, 
#' and conduct sex assignment. Molecular Ecology Resources.
#' }
#' @examples
#' filtered.gl <- gl.filter.excess.het(x = LBP, Yates = TRUE)
#' # Use below function to output information of the loci with Yates's continuity correction specified 
#' filtered.table <- gl.report.excess.het(x = LBP, Yates = TRUE)
#' @seealso \code{\link{gl.filter.callrate}}
#' @importFrom stats aggregate
#' @export
#' @return Returns unaltered genlight object

gl.filter.excess.het <- function(x,
                                 Yates = FALSE,
                                 mono.rm = FALSE,
                                 recalc = FALSE,
                                 verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # DO THE JOB
  
  if (!Yates) {
    cc = 0
  } else {
    cc = 0.5
  }
  
  # Start calculating stats from observed data
  # Plot
  if (verbose >= 2) {
    cat(report("  Calculating loci statistics from observed data"))
  }
  
  gen <- as.data.frame(t(as.matrix(x)))
  n0 <- rowSums(gen == 0, na.rm = TRUE)
  n1 <- rowSums(gen == 1, na.rm = TRUE)
  n2 <- rowSums(gen == 2, na.rm = TRUE)
  
  fhe <- n1 / (n0 + n1 + n2)
  
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
    
  # OUTPUT FILTERED GENLIGHT
  # Remove highly-het loci from new filtered gl
   x2 <- x[,!(x$loc.names %in% table.filter$loci)]
  
  # Remove highly-het loci from new gl loci metadata
   x2@other$loc.metrics <-
    x@other$loc.metrics[!(x$loc.names %in% table.filter$loci),]
  
  # Output
  if (verbose >= 2) {
    cat(report(
      "  Removed ",
      x@n.loc - x2@n.loc,
      " excessively-heterozygous loci."
    ))
  }
  
  # Monomorphic loci may have been created
  x2@other$loc.metrics.flags$monomorphs <- FALSE
  
  # Monomorphic loci may have been created ------
  if (nrow(table.filter) != 0 & mono.rm==TRUE) {
    if (verbose >= 2) {
      cat(report("  Deleting monomorphic loc\n"))
    }
    x2 <- gl.filter.monomorphs(x2, verbose = 0)
  }
  
  # Check monomorphs have been removed
  if (nrow(table.filter) != 0 & x2@other$loc.metrics.flags$monomorphs == FALSE) {
    if (verbose >= 2) {
      cat(warn(
        "  Warning: Resultant dataset may contain monomorphic loci\n"
      ))
    }
  }
  
  # Recalculate statistics ----------
  if (nrow(table.filter) != 0 & recalc) {
    x2 <- gl.recalc.metrics(x2, verbose = 0)
    if (verbose >= 2) {
      cat(report("  Recalculating locus metrics\n"))
    }
  } else if (nrow(table.filter) != 0 & recalc==F) {
    if (verbose >= 2) {
      cat(warn("  Locus metrics not recalculated\n"))
      x2 <- utils.reset.flags(x2, verbose = 0)
    }
  }
  
  # ADD TO HISTORY
  
  nh <- length(x2@other$history)
  x2@other$history[[nh + 1]] <- match.call()
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(invisible(x2))
  
}
