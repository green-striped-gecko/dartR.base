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
#' entities and individuals as attributes.
#' 
# @export
#' @return a transposed genlight object

utils.transpose <- function(x,
                            parallel = FALSE) {
  hold <- x
  # Store the input object x in a temporary variable hold
  
  x@gen <- matrix2gen(t(as.matrix(x)), parallel = parallel)
  # Transpose the genotype matrix of x using the t() function and convert it to a genind object using matrix2gen()
  # The transposed genotype matrix is assigned to the gen slot of the x object
  
  x@n.loc <- nInd(hold)
  # Set the number of loci in x to the number of individuals in the original object x
  
  indNames(x) <- locNames(hold)
  locNames(x) <- indNames(hold)
  # Swap the individual names and locus names in x with the corresponding names from the original object
  
  # This is just a dummy vector to comply with the attributes of a genlight object
  alleles(x) <- paste(rep("A", nInd(hold)), rep("A", nInd(hold)), sep = "/")
  # Set the alleles slot of x to a dummy vector consisting of "A/A" repeated for the number of individuals
  
  ploidy(x) <- unique(ploidy(hold))
  # Set the ploidy of x to the unique ploidy value in the original object
  
  pop(x) <- rep("NA", nLoc(hold))
  # Set the population slot of x to "NA" repeated for the number of loci
  
  x@other$loc.metrics <- hold@other$ind.metrics
  x@other$ind.metrics <- hold@other$loc.metrics
  # Swap the values of the loc.metrics and ind.metrics slots in the other slot of x with the corresponding values from the original object
  
  return(x)
  # Return the modified x object
}