#' @name gl.filter.replicates
#' @title Remove replicated individuals 
#' @description
#' Remove replicated individuals 
#' @param x Name of the genlight object containing the SNP data [required].
#' @param replicates.report Output of function gl.report.replicates [required].
#' @param loc_threshold Minimum number of loci required to assess that two
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
#' gl.report.replicates. Note that thresholds can only be tightened relative
#' to those used in gl.report.replicates -- the report table contains only
#' the pairs that passed the original thresholds. For each retained pair the
#' member with the higher missing-data rate is removed; ties remove the
#' alphabetically later individual.
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \donttest{
#' t1 <- platypus.gl[1:50,]
#' res_rep <- gl.report.replicates(t1, loc_threshold = 500, 
#' perc_geno = 0.85)
#' t2 <- gl.filter.replicates(t1, replicates.report = res_rep, perc_geno = 0.85)
#' }
#' @family matched filter
#' @return A reduced dartR genlight object, returned invisibly
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

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # Validate the report input
  if (!is.list(replicates.report) ||
      !("table.rep" %in% names(replicates.report)) ||
      !is.data.frame(replicates.report$table.rep)) {
    stop(error(
      "Fatal Error: replicates.report must be the list returned by gl.report.replicates\n"
    ))
  }

  RR <- data.table::as.data.table(replicates.report$table.rep)

  col_same <- RR[nloc > loc_threshold & perc > perc_geno][order(-perc)]

  if (nrow(col_same) > 0) {
    # Canonicalise pair orientation (reports generated before the pair
    # table was deduplicated carry both orderings of each pair) and
    # reduce to one row per unordered pair
    swap <- col_same$ind1 > col_same$ind2
    if (any(swap)) {
      tmp.ind <- col_same$ind1[swap]
      tmp.miss <- col_same$ind1_miss[swap]
      col_same$ind1[swap] <- col_same$ind2[swap]
      col_same$ind2[swap] <- tmp.ind
      col_same$ind1_miss[swap] <- col_same$ind2_miss[swap]
      col_same$ind2_miss[swap] <- tmp.miss
    }
    col_same <- col_same[!duplicated(paste(ind1, ind2))]

    # Decide which replicate to drop: the one with the higher missing
    # proportion; ties drop ind2 (the alphabetically later individual)
    col_same[, ind_to_drop := ifelse(ind1_miss > ind2_miss, ind1, ind2)]
  }
  drop_list <- unique(col_same$ind_to_drop)

  if (length(drop_list) == 0) {
    if (verbose >= 1) {
      cat(report(
        "  No replicate pairs pass the thresholds; no individuals removed\n"
      ))
      cat(report("Completed:", funname, "\n"))
    }
    return(invisible(x))
  }

  if (verbose >= 2) {
    cat(report(
      "  Removing replicated individuals:",
      paste(drop_list, collapse = ", "),
      "\n"
    ))
  }

  hold.history <- x@other$history
  xx <- gl.drop.ind(x,
                    ind.list = drop_list,
                    recalc = recalc,
                    mono.rm = mono.rm,
                    verbose = 0)

  # ADD TO HISTORY (a single entry for this call)
  xx@other$history <- hold.history
  nh <- length(xx@other$history)
  xx@other$history[[nh + 1]] <- match.call()

  # Final verbose message if requested
  if (verbose >= 1) cat(report("Completed:", funname, "\n"))

  return(invisible(xx))

}
