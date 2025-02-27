#' @name utils.collapse.matrix
#' @title Collapses a distance matrix calculated for individuals to a distance matrix for populations defined in a dartR
#' genlight object
#' @family utilities

#' @description 
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.

#' @param D Name of the matrix containing the distances between individuals [required].
#' @param x Name of the genlight object containing the genotypes [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log; 3, progress and results summary; 5, full report [default 2].
#'  
#'  @details
#' This script takes a matrix of distances calculated between individuals and collapses it by averaging to a matrix of distances between populations.
#' The script gl.dist.ind has a lot of options for distances for presence absence data, but gl.dist.pop does not. This script allows efficient and
#' consistent transfer of this capability to gl.dist.pop.

#' @author Author: Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @export  
#' @return An object of class 'dist' or 'matrix' giving distances between individuals

# Examples for testing
# D.ind <- gl.dist.ind(testset.gl)
# as.matrix(D.ind)[1:7,1:7]
# D.pop <- utils.collapse.matrix(D=D.ind,x=testset.gl,verbose=3)
# as.matrix(D.pop)[1:3,1:3]

utils.collapse.matrix <- function(D,
                                  x,
                                  verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype1 <-
    utils.check.datatype(x, accept = c("SNP","SilicoDArT"),verbose = verbose)
  
  if (!is(x, "dartR")) {
    class(x) <- "dartR"  
    if (verbose>2) {
      cat(warn("Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"))
      cat(warn("                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR, 
                 please use function dartR2gl\n"))
    }
  }
  
  datatype2 <-
    utils.check.datatype(D, accept = c("dist","matrix"), verbose = verbose)
  
  # SCRIPT SPECIFIC ERROR CHECKING
  
  if(datatype2=="dist"){
    mat <- as.matrix(D)
  } else {
    mat <- D
  }
  
  # Check if matrix is square
  if (!is.matrix(mat) || nrow(mat) != ncol(mat)) {
    cat(error("Fatal Error: Input must be a square distance matrix or object of class 'dist'\n"))
    stop()
  }
  
  # Ensure row and column names in the matrix match individual names
  if (!all(rownames(mat) %in% indNames(x)) || !all(colnames(mat) %in% indNames(x))) {
    cat(error("Fatal Error: Matrix row/column names do not match individual names in genlight object\n"))
    stop()
  }
  
  # DO THE JOB
  
  # Create an empty matrix for population-level distances
  pop.mat <- matrix(0, nrow = nPop(x), ncol = nPop(x), dimnames = list(popNames(x), popNames(x)))
  
  # Compute average pairwise distances between populations
  if(verbose >= 3){cat(report("  Computing mean distances between populations\n"))}
  for (i in seq_along(popNames(x))) {
    for (j in seq_along(popNames(x))) {
      # Get individuals belonging to each population
      inds_i <- indNames(x)[pop(x) == popNames(x)[i]]
      inds_j <- indNames(x)[pop(x) == popNames(x)[j]]
      
      # Extract the relevant submatrix
      submat <- mat[inds_i, inds_j, drop = FALSE]
      
      # Compute mean distance, handling cases where populations are the same
      if (length(submat) > 0) {
        pop.mat[i, j] <- mean(submat, na.rm = TRUE)  # Average distance
      }
    }
  }
  
  if(class(D)=='dist'){
    pop.mat <- as.dist(pop.mat)
    if(verbose >= 3){cat(report("  Returning object of class 'dist'\n"))}
  }  else {
    if(verbose >= 3){cat(report("  Returning object of class 'matrix'\n"))}
  }
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(pop.mat)
}

