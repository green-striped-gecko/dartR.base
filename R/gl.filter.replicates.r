#' @name gl.filter.replicates
#' @title Remove replicated individuals 
#' @description
#' Remove replicated individuals 
#' @param x Name of the genlight object containing the SNP data [required].
#' @param replicates.report Output of functionv gl.report.replicates [required].
#' @param loc_threshold Minimum number of loci required to asses that two 
#' individuals are replicates [default 100].
#' @param perc_geno Minimum percentage of genotypes in which two individuals 
#' should be the same [default 0.95]. 
#' @param recalc If TRUE, recalculate the locus metadata statistics 
#' [default FALSE].
#' @param mono.rm If TRUE, remove monomorphic and all NA loci [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' Remove replicated individuals using as input the report from function 
#' gl.report.replicates. The user can choose new thresholds for the minimum 
#' number of loci required to asses that two individuals are replicates 
#' (loc_threshold) and the minimum percentage of genotypes in which two 
#' individuals should be the same (perc_geno) from those thresholds use in 
#' gl.report.replicates. 
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' t1 <- platypus.gl
#' res_rep <- gl.report.replicates(t1, loc_threshold = 500, 
#' perc_geno = 0.85)
#' t2 <- gl.filter.replicates(t1, replicates.report = res_rep, perc_geno = 0.85)
#' @family matched filter
#' @return A reduced dartR genlight object
#' @export

gl.filter.replicates <- function(x,
                                 replicates.report, 
                                 loc_threshold = 100,    
                                 perc_geno = 0.95,
                                 recalc = FALSE,
                                 mono.rm = FALSE,
                                 verbose = NULL){
  
  ind1 <- ind1_miss <- ind2 <- ind2_miss <- ind_to_drop <- nloc <- perc <- NULL
  
  # Determine verbosity level (internal helper)
  verbose  <- gl.check.verbosity(verbose)
  
  # Record function name for logging
  funname  <- match.call()[[1]]
  
  # Flag the start of this function call (internal helper)
  utils.flag.start(func = funname, build = "Jody", verbose = verbose)
  
  RR <- replicates.report$table.rep
  
  col_same <- RR[nloc > loc_threshold & perc > perc_geno][order(-perc)]
  
  # Decide which replicate to drop: the one with higher missing proportion
  col_same[, ind_to_drop := ifelse(ind1_miss > ind2_miss, ind1, ind2)]
  drop_list <- unique(col_same$ind_to_drop)
  
  xx <- gl.drop.ind(x, 
                    ind.list = drop_list,
                    recalc = recalc,
                    mono.rm = mono.rm,
                    verbose = verbose)
  

  # Final verbose message if requested
  if (verbose >= 1) cat(report("Completed:", funname, "\n"))
  
  return(xx)

}
