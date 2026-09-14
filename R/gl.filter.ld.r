#' @name gl.filter.ld
#' @title Filters loci based on linkage disequilibrium (LD)
#' @family matched filter
#' 
#' @description
#' This function uses the statistic set in the parameter \code{stat.keep} from
#' function \code{\link{gl.report.ld.map}} to choose the SNP to keep when two
#' SNPs are in LD. When a SNP is selected to be filtered out in a pairwise
#' comparison, the function stores its name in a list; a SNP already in the
#' list is not listed again.
#'
#' @details
#' Within each population, the pairs reported by
#' \code{\link{gl.report.ld.map}} with \code{ld.stat} at or above
#' \code{threshold} are processed in the order they appear in the report. For
#' each pair, the locus with the lower \code{stat.keep} value is marked for
#' removal. Note that comparisons are strictly pairwise and sequential: a
#' locus can be marked for removal against a partner that was itself marked
#' by an earlier comparison. A locus is removed from the dataset when it is
#' marked in at least \code{pop.limit} populations.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param ld.report Output from function \code{\link{gl.report.ld.map}}
#' [required].
#' @param threshold Threshold value at or above which loci will be removed
#' [default 0.2].
#' @param pop.limit Minimum number of populations in which LD should be more
#' than the threshold for a locus to be filtered out [default half of the
#' populations represented in ld.report, i.e. ceiling(n/2)].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @author Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' test <- gl.filter.callrate(platypus.gl, threshold = 1)
#' test <- gl.filter.monomorphs(test)
#' res <- gl.report.ld.map(test)
#' res_2 <- gl.filter.ld(x = test, ld.report = res)
#' }
#' @seealso \code{\link{gl.report.ld.map}}
#'
#' @export
#' @return The reduced genlight object (returned invisibly).

gl.filter.ld <- function(x,
                         ld.report,
                         threshold = 0.2,
                         pop.limit = NULL,
                         verbose = NULL) {

  x_hold <- x

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  # ld.report must be the data frame produced by gl.report.ld.map
  req.cols <- c("pop", "ld.stat", "locus_a.snp.name", "locus_b.snp.name",
                "locus_a.stat.keep", "locus_b.stat.keep")
  if (!is.data.frame(ld.report) ||
      !all(req.cols %in% colnames(ld.report))) {
    stop(error(
      "Fatal Error: ld.report must be the data frame returned by",
      " gl.report.ld.map; required columns:",
      paste(req.cols, collapse = ", "), "\n"
    ))
  }

  # default pop.limit: half of the populations represented in the report
  # (evaluated explicitly here; populations skipped by gl.report.ld.map's
  # ind.limit are not counted)
  if (is.null(pop.limit)) {
    pop.limit <- ceiling(length(unique(ld.report$pop)) / 2)
  }

  # Check monomorphs have been removed up to date
  if (!isTRUE(x@other$loc.metrics.flags$monomorphs)) {
    if (verbose >= 2) {
      cat(
        warn(
          "  Warning: Data may include monomorphic loci in call rate
          calculations for filtering\n"
        )
      )
    }
  }

  x <- gl.keep.pop(x,pop.list = as.character(unique(ld.report$pop)),verbose = 0)

  ld_tmp <- ld.report[ld.report$ld.stat >= threshold, ]
  if(nrow(ld_tmp) == 0){
    if (verbose >= 1) {
      cat(report(paste(
        " No pair of loci were found to be in LD using a threshold of",
        threshold,"\n","Returning an unaltered genlight object\n ")))
    }
    # FLAG SCRIPT END

    if (verbose > 0) {
      cat(report("Completed:", funname, "\n"))
    }

    return(invisible(x_hold))
  }else{
  ld_tmp$test_stat <- ld_tmp$locus_a.stat.keep >= ld_tmp$locus_b.stat.keep
  ld_tmp$pop <- as.factor(ld_tmp$pop)
  ld_tmp_pop <- split(ld_tmp, f = ld_tmp$pop)
  
  loci_list <- vector(mode = "list", length = length(ld_tmp_pop))
  
  for (i in 1:length(ld_tmp_pop)) {
    ld_pop <- ld_tmp_pop[[i]]
    for (y in 1:nrow(ld_pop)) {
      if (ld_pop[y, "test_stat"] == TRUE) {
        loci_tmp <- ld_pop[y, "locus_b.snp.name"]
      } else{
        loci_tmp <- ld_pop[y, "locus_a.snp.name"]
      }
      
      if (loci_tmp %in% loci_list[[i]]) {
        next
      } else{
        loci_list[[i]] <- c(loci_list[[i]], loci_tmp)
      }
      
    }
  }
  
  loci_list_res <- Reduce("c", loci_list)
  loci_names_tmp <- names(table(loci_list_res))
  loci_names <- loci_names_tmp[table(loci_list_res) >= pop.limit]
  
  x2 <- gl.drop.loc(x_hold, loc.list =  loci_names, verbose = 0)

  # the internal gl.drop.loc call appends its own history entry exposing
  # internal variable names; reset so this function adds the single entry
  x2@other$history <- x_hold@other$history

  # REPORT A SUMMARY
  if (verbose >= 2) {
    cat("  Summary of filtered dataset\n")
    cat(paste("    LD for loci >", threshold, "\n"))
    cat(paste("    Original No. of loci :", nLoc(x), "\n"))
    cat(paste("    No. of loci retained:", nLoc(x2), "\n"))
    cat(paste("    No. of populations: ", nPop(x2), "\n"))
  }
  
  # ADD TO HISTORY
  
  nh <- length(x2@other$history)
  x2@other$history[[nh + 1]] <- match.call()
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(invisible(x2))
}
  
}
