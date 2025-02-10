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
#' @param by.pop If FALSE, ignore population settings when subsampling; if TRUE, subsample
#' each population to n individuals [default TRUE]. 
#' @param error.check If TRUE, will undertake error checks on input paramaters [default TRUE]
#' @param mono.rm If TRUE and error.check is TRUE, monomorphic loci arising from the deletion of individuals
#' will be filtered from the resultant genlight object [default FALSE]
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
                  mono.rm=FALSE,
                  verbose = NULL) {
  
  if(error.check==TRUE){
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose=verbose)
    
    if (!is(x, "dartR")) {
      class(x) <- "dartR"  
      if (verbose>2) {
        cat(warn("Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"))
        cat(warn("                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR, 
                 please use function dartR2gl\n"))
      }
    }
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if(n < 1){
      cat(error("Fatal Error: Number of individuals must be a positive integer\n"))
      stop()
    }
    
    # Set the limit for n if Replace=FALSE
    n.limit <- min(table(pop(x)))
    
    # Set the default value of n if not specified
    if (is.null(n)){
      if (by.pop){
        n <- n.limit
      } else {
        n <- nInd(x)
      }
    }
    
    if (by.pop==TRUE){
      if ( n > n.limit){
        if (verbose >= 1) {cat(warn("  Warning: Specified subsample size larger than smallest population size\n"))}
        if(replace==FALSE){
          if (verbose >= 1) {cat(warn("    Sampling without replacement, n set to",n.limit,"\n"))}
          n <- n.limit
        } else {
          if (verbose >= 1) {cat(warn("    Some populations will be upsampled by replacement\n"))}
        }
      }
    }
    if(by.pop==FALSE){
      if ( n > nLoc(x) ){
        if (verbose >= 1) {cat(warn("  Warning: Specified subsample size larger than total number of individuals\n"))}
        if(replace==FALSE){
          if (verbose >= 1) {cat(warn("    Sampling without replacement, n set to",nLoc(x),"\n"))}
          n <- nLoc(x)
        } else {
          if (verbose >= 1) {cat(warn("    Some populations will be upsampled by replacement\n"))}
        }
        n <- nLoc(x)
      }
    }
  }
  # DO THE JOB

  # Internal function
  subsample <- function(x,n,replace){
    if(n <= nInd(x)){
      idx <- sample(seq_len(nInd(x)), size = n, replace = replace)
      xx <- x[idx,] # Requires dartR genlight object to ensure loc.metrics subsetted also
    } 
    if(n > nInd(x)){
      if(replace==TRUE){
        # Subsample with replacement to the existing number of individuals
        idx <- sample(seq_len(nInd(x)), size = nInd(x), replace = TRUE)
        xx <- x[idx,] # Requires dartR genlight object to ensure loc.metrics subsetted also
        # If the requested sample sizes are greater than nInd(x)
        # then if n is more than double nInd(x)
        if(n/nInd(x)>=2){
          # Add increments of nInd(x) individuals to the new genlight object
          for (i in 1:trunc(n/nInd(x))-1){
            idx <- sample(seq_len(nInd(x)), size = nInd(x), replace = TRUE)
            tmp <- x[idx,] # Requires dartR genlight object to ensure loc.metrics subsetted also
            xx <- gl.join(xx,tmp,method="end2end",verbose=0)
          }
        }
        # Deal with the remainder if any
        if(n%%nInd(x) > 0){
          idx <- sample(seq_len(nInd(x)), size = (n %% nInd(x)), replace = TRUE)
          tmp <- x[idx,]
          xx <- gl.join(xx,tmp,method="end2end",verbose=0)
        }
      } else {
        cat(error("Fatal Error: Cannot upsample a genlight object without replacement\n"))
        stop()
      }
    }
    return(xx)
  }

  
  if(by.pop==FALSE){
    xx <- subsample(x=x,n=n,replace=replace)
   }
  
   if(by.pop == TRUE){
    new.list <- list()  # Store subsampled genlight objects
    
    # Iterate over population names to subset and subsample
    pop_names <- popNames(x)  # Get population names
    for(i in seq_along(pop_names)){
      gl <- gl.keep.pop(x, pop.list = pop_names[i], verbose = 0)
      new.list[[i]] <- subsample(x = gl, n = n, replace = replace)
    }
    
    # Ensure the list is not empty
    if(length(new.list) > 1){
      xx <- new.list[[1]]  # Start with the first population's subset
      for(i in 2:length(new.list)){
        xx <- gl.join(xx, new.list[[i]], method="end2end",verbose=0)
      }
    } else {
      xx <- new.list[[1]]  # If only one population, return as is
    }
  }

    if(error.check==TRUE){
      # FILTER MONOMORPHS
      if(mono.rm==TRUE){
        x <- gl.filter.monomorphs(x,verbose=verbose)
        
        # ADD TO HISTORY
        nh <- length(x2@other$history)
        x2@other$history[[nh + 1]] <- match.call()
      }
      
      # FLAG SCRIPT END ---------------
      
      if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
      }
    }
  
  return(xx)
}
    
