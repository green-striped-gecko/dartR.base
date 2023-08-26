#' @name gl.subsample.loc
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
#' @param error.check If TRUE, will undertake error checks on input paramaters [default TRUE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#'  
#' @details Retain a subset of loci at random, with or without replacement.
#' Parameter n must be less than or equal to nLoc(x). 
#' 
#' #' Set error.check = FALSE for speedy execution in simulations
#' 
#' @author Custodian: Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})

#' @examples 
#' gl2 <- gl.subsample.loc(testset.gl, n=50, replace=TRUE, verbose=3)

#' @export 
#' @return Returns the subsampled genlight object

gl.subsample.loc <- function(x,
                              n,
                              replace=TRUE,
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
    datatype <- utils.check.datatype(x, verbose = verbose)

    # FUNCTION SPECIFIC ERROR CHECKING

    if (n <= 0 | n > nLoc(x)) {
      if (verbose >= 1) {cat(warn("Subsample size must be in the range 1 to",nLoc(x),"\n"))}
      if (verbose >= 1) {cat(warn("  Set to",nLoc(x),"\n"))}
      n <- nLoc(x)
    }

    if (verbose >= 2) {
      if (replace){
        if (verbose >= 2) {cat(report("  Subsampling",n,"loci at random from a",datatype,"object with replacement\n"))}
      } else {
        if (verbose >= 2) {cat(report("  Subsampling",n,"loci at random from a",datatype,"object without replacement\n"))}
      }
    }
  }
    
  # DO THE JOB
    
  # Subsample the genlight object
  # generate a random index value, with or without replacement
  nums <- sample(1:nLoc(x), size = n, replace = replace)
  # subsample the data
  x2 <- x[,nums]
  # subsample the locus metrics [necessary because of replacment possibility]
  x2@other$loc.metrics <- x@other$loc.metrics[nums,]

  # if(error.check){
  # # ADD TO HISTORY
  #   nh <- length(x2@other$history)
  #   x2@other$history[[nh + 1]] <- match.call()
  #   
  # # FLAG SCRIPT END ---------------
  #   
  #   if (verbose >= 1) {
  #     cat(report("Completed:", funname, "\n"))
  #   }
  # } 
    
  return(x2)
}
    
