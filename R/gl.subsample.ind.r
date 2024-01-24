#' @name gl.subsample.ind
#' @title Subsample individuals from a genlight object
#' @family data manipulation
#'
#' @description 
#' A function to subsample individuals at random in a genlight object
#' with and without replacement.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param n Number of individuals to include in the subsample [default NULL]
#' @param replace If TRUE, sampling is with replacement [default TRUE]
#' @param by.pop If FALSE, ignore population settings [default TRUE]. 
#' @param error.check If TRUE, will undertake error checks on input paramaters [default TRUE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#'  
#' @details Retain a subset of individuals at random, with or without replacement.
#' If subsampling globally, n must be less than or equal to nInd(x). If subsampling
#' by population, then n must be less than the minimum sample size for any population.
#' 
#' Set error.check = FALSE for speedy execution in simulations
#' 
#' @author Custodian: Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})

#' @examples 
#' gl <- gl.subsample.ind(testset.gl, n=30, by.pop=FALSE, replace=TRUE)
#' gl <- gl.subsample.ind(platypus.gl, n=10, by.pop=TRUE, replace=TRUE)

#' @export 
#' @return Returns the subsampled genlight object

gl.subsample.ind <- function(x,
                  n = NULL,
                  replace = TRUE,
                  by.pop = TRUE,
                  error.check = TRUE,
                  verbose = NULL) {
  
  if(error.check){
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose=verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    n.limit <- min(table(pop(x)))
    if (is.null(n)){
      if (by.pop){
        n <- n.limit
      } else {
        n <- nInd(x)
      }
    }
    
    if (by.pop){
      if ( n > n.limit){
        if (verbose >= 1) {cat(warn("  Warning: Specified subsample size larger than smallest population size\n"))}
        if (verbose >= 1) {cat(warn("    n set to",n.limit,"\n"))}
        n <- n.limit
      }
    } else {
      if ( n > nLoc(x) ){
        if (verbose >= 1) {cat(warn("  Warning: Specified subsample size larger than total number of individuals\n"))}
        if (verbose >= 1) {cat(warn("    n set to",nLoc(x),"\n"))}
        n <- nLoc(x)
      }
    }
  }
  # DO THE JOB

  if (!by.pop){
    # Generate a random set of n numbers
      nums <- sample(1:nInd(x), size = n, replace = replace)
      ind.list <- indNames(x)[nums]
    # Subsample the genlight object
      x2 <- gl.keep.ind(x,ind.list=ind.list,verbose=0)
    #x2@other <- x@other[nums,]
  } else {
    popcount <- 1
    for (popn in popNames(x)){
      tmp <- gl.keep.pop(x,pop.list=popn,verbose=0)
      # Generate a random set of n numbers
      nums <- sample(1:nInd(tmp), size = n, replace = replace)
      ind.list <- indNames(tmp)[nums]
      tmp <- gl.keep.ind(tmp,ind.list=ind.list,verbose=0)
      #tmp <- tmp[nums,]
      if(popcount == 1){
        x2 <- tmp
      } else {
        hold <- x2@other$ind.metrics
        x2 <- rbind(x2,tmp)
        x2@other$ind.metrics <- rbind(hold,tmp@other$ind.metrics)
      }
      popcount <- popcount + 1
    }
  }
  x2@other$loc.metrics <- x@other$loc.metrics

  if(error.check){
    # ADD TO HISTORY
    nh <- length(x2@other$history)
    x2@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END ---------------
    
    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }
  }
  
  return(x2)
}
    
