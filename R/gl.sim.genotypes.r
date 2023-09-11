#' @name gl.sim.genotypes
#' @title Generate random genotypes 
#' @family data manipulation
#' 
#' @description
#' Generate random genotypes for a single population drawing upon the allele
#' frequencies from that population.
#' 
#' @param x Name of the genlight object [required].
#' @param n.ind Number of individuals to be simulated (should be less than the number of loci) [default 200]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity]
#' 
#' @author Custodian: Arthur Georges (Post to
#' \url{https://groups.google.com/d/forum/dartr}) 
#' 
#' @export 
#' @return Returns a genlight object with the simulated genotypes

gl.sim.genotypes <- function(x,
                             n.ind=200,
                             #error.check = TRUE,
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
  
  
  # Extract the allele frequencies from the given genlight object
  df <- gl.allele.freq(x,by="loc",verbose=0)
  n.loc <- nrow(df)
  
  # Warning on more individuals than loci
  if(n.ind > n.loc){
    cat(warn("  Warning: the number of individuals (entities) exceeds the number of loci (attributes) which will cause issues for analyses like PCA\n"))
    cat(warn("    Setting n.ind to",n.loc,"\n"))
  }

  # Create an array to hold the new genotypes
  v <- array(NA,dim=c(n.loc,n.ind))
  
  # Populate the array  
  for (i in 1:n.loc){
    # Generate a vector of individuals, sampling locus i using the allele frequencies for that locus
    v1 <- sample(c(0, 1), size = n.ind, replace = TRUE, prob = c((1-df$frequency[i]),df$frequency[i]))
    # Sample again
    v2 <- sample(c(0, 1), size = n.ind, replace = TRUE, prob = c((1-df$frequency[i]),df$frequency[i]))
    # Combine the two random haplotypes to form a genotype
    v[i,] <- v1 + v2
  }
  
  # Convert the array v to a new genlight object  
  gl <- new(
    "genlight",
    gen = t(v),
    ind.names = paste0("Ind_", 1:n.ind),
    loc.names = as.character(df[,1]),
    ploidy = rep(2, n.loc)
  )
  
  # Enforce its compliance with dartR
  gl <- gl.compliance.check(gl,verbose=0)
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(gl)
  
}
