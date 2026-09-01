#' @name utils.transpose
#' @title An internal utility function to transpose a genlight object.
#' @family utilities
#' 
#' @description 
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.

# build = "v.2023.2"

#' @param x name of the genlight object
#' @param parallel if TRUE, use parallel processing capability
#' 
#' @details
#' This is a function to transpose a genlight object, that is, to set loci as
#' entities and individuals as attributes. Depends on matrix2gen
#' (utils.impute.R). Assumes uniform ploidy.
#' 
#' @keywords internal
# @export
#' @return a transposed genlight object

utils.transpose <- function(x,
                            parallel = FALSE) {
  hold <- x
  
  x@gen <- matrix2gen(t(as.matrix(x)), parallel = parallel)
  
  x@n.loc <- nInd(hold)
  
  indNames(x) <- locNames(hold)
  locNames(x) <- indNames(hold)
  
  # This is just a dummy vector to comply with the attributes of a genlight object
  alleles(x) <- paste(rep("A", nInd(hold)), rep("A", nInd(hold)), sep = "/")
  
  ploidy(x) <- unique(ploidy(hold))
  
  pop(x) <- rep("NA", nLoc(hold))
  
  x@other$loc.metrics <- hold@other$ind.metrics
  x@other$ind.metrics <- hold@other$loc.metrics
  
  return(x)
}