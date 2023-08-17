#' @name gl.subsample.loc.v2
#' @title Subsample loci from a genlight object
#' @family data manipulation
#' 
#' @description 
#' A function to subsample loci at random in a genlight object
#' with and without replacement.

#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param n Number of loci to include in the subsample [default NULL]
#' @param replace If TRUE, sampling is with replacement [default TRUE]
#' @param mono.rm If TRUE, monomorphic loci are excluded [default TRUE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#'  
#' @details Retain a subset of iloci at random, with or without replacement.
#' Parameter n must be less than or equal to nLoc(x). 
#' 
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})

#' @examples 
#' gl <- gl.subsample.loc.v2(testset.gl, n=100, replace=TRUE)

#' @export 
#' @return Returns the subsampled genlight object

gl.subsample.loc.v2 <- function(x,
                              n,
                              replace=TRUE,
                              mono.rm = TRUE,
                              verbose = NULL) {

    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    # Check for monomorphs
    if (!mono.rm){
      if (x@other$loc.metrics.flags$monomorphs == FALSE) {
        if (verbose >= 2) {
          cat(warn("  Warning: Dataset may contain monomorphic loci\n"))
        }
      }
    } else {
      x <- gl.filter.monomorphs(x, verbose = 0)
      if (verbose >= 2) {cat(report("  Removing monomorphic loci\n"))}
    }
    
    if (n <= 0 | n > nLoc(x)) {
      if (verbose >= 1) {cat(warn("Subsample size must be in the range 1 to",nLoc(x),"\n"))}
      if (verbose >= 1) {cat(warn("  Set to",nLoc(x),"\n"))}
        n <- nLoc(x)
    }
    
    # DO THE JOB
    
    if (verbose >= 2) {
      if (replace){
        if (verbose >= 2) {cat(report("  Subsampling",n,"loci at random from a",datatype,"object with replacement\n"))}
        nums <- sample(1:nLoc(x), size = n, replace = TRUE)
      } else {
        if (verbose >= 2) {cat(report("  Subsampling",n,"loci at random from a",datatype,"object without replacement\n"))}
        nums <- sample(1:nLoc(x), size = n, replace = FALSE)
      }
    }
    # Subsample the genlight object
    x2 <- x[,nums]
    
    # ADD TO HISTORY
    nh <- length(x2@other$history)
    x2@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END ---------------
    
    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }
    
    return(x2)
}
    
