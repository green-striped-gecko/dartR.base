#' @name gl.filter.replicates.KING
#' @title Removes duplicate samples identified by KING-robust kinship
#' @family matched filter
#'
#' @description
#' Removes the duplicate samples identified by
#' \code{\link{gl.report.replicates.KING}}, keeping one sample of each group
#' of duplicates.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param replicates.report Output of \code{\link{gl.report.replicates.KING}}
#' [required].
#' @param king_threshold KING-robust kinship above which a pair is a
#' duplicate [default 0.354].
#' @param het_threshold Heterozygote share above which a pair is a duplicate;
#' used only when the report was computed from called genotypes
#' [default 0.80].
#' @param loc_threshold Minimum number of loci two samples must share for
#' their pair to be assessed [default 100].
#' @param recalc If TRUE, recalculate the locus metadata statistics
#' [default FALSE].
#' @param mono.rm If TRUE, remove monomorphic and all NA loci [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' The pairs are flagged again from \code{replicates.report$table.all} with
#' the thresholds given here: a pair is a duplicate when it shares at least
#' \code{loc_threshold} loci and its KING-robust kinship exceeds
#' \code{king_threshold} or, when the report was computed from called
#' genotypes (\code{replicates.report$input}), its heterozygote share exceeds
#' \code{het_threshold}, as in the report. Because table.all holds every pair, the thresholds
#' can be loosened as well as tightened, except \code{loc_threshold}: pairs
#' below the report's \code{loc_threshold} have no statistics. Flagged pairs
#' are linked into groups and, in each group, the sample with the lowest
#' proportion of missing genotypes (\code{replicates.report$ind.miss}) is
#' kept (ties keep the name that sorts first) and the others are removed.
#' With the report's thresholds the filter removes exactly
#' \code{replicates.report$ind.list.drop}.
#'
#' Run \code{\link{gl.filter.contamination}} before this function. Keeping
#' the sample with the fewest missing genotypes favours a contaminated
#' sample, as mixed DNA gives more calls: in a ddRAD study of collared
#' peccaries the pair CP47-CP57 keeps CP57, a mixture of about 4/5 CP47 and
#' 1/5 CP46.
#'
#' Pairs inside a group that are not flagged themselves are listed with a
#' warning: such a group can link two animals through a sample that is a
#' mixture of them, and removing all but one sample would then remove a
#' different animal. See \code{\link{gl.report.replicates.KING}}.
#'
#' @return The genlight object with the duplicate samples removed, returned
#' invisibly.
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#' t1 <- platypus.gl[1:50, ]
#' res <- gl.report.replicates.KING(t1, plot.out = FALSE)
#' t2 <- gl.filter.replicates.KING(t1, replicates.report = res)
#' }
#'
#' @seealso \code{\link{gl.report.replicates.KING}},
#' \code{\link{gl.filter.replicates}}
#' @export

gl.filter.replicates.KING <- function(x,
                                      replicates.report,
                                      king_threshold = 0.354,
                                      het_threshold = 0.80,
                                      loc_threshold = 100,
                                      recalc = FALSE,
                                      mono.rm = FALSE,
                                      verbose = NULL) {

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, build = "Jody", verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (!is.list(replicates.report) ||
      !all(c("table.all", "ind.miss", "input") %in% names(replicates.report)) ||
      !is.data.frame(replicates.report$table.all) ||
      !isTRUE(replicates.report$input %in% c("likelihoods", "counts", "calls"))) {
    stop(error("  replicates.report must be the list returned by gl.report.replicates.KING\n"))
  }
  number <- function(v) is.numeric(v) && length(v) == 1 && !is.na(v)
  if (!number(king_threshold) || !number(het_threshold) ||
      !number(loc_threshold)) {
    stop(error("  king_threshold, het_threshold and loc_threshold must be numbers\n"))
  }
  rep.loc <- replicates.report$parameters$loc_threshold
  if (!is.null(rep.loc) && loc_threshold < rep.loc && verbose >= 1) {
    cat(warn("  Warning: the report has no statistics for pairs sharing fewer than",
             rep.loc, "loci, so they cannot be flagged\n"))
  }

  # The heterozygote share flags pairs only from called genotypes, as in
  # the report
  use.het <- replicates.report$input == "calls"
  if (verbose >= 2) {
    if (use.het) {
      cat(report("  Flagging pairs with KING above", king_threshold,
                 "or a heterozygote share above", het_threshold, "\n"))
    } else {
      cat(report("  Flagging pairs with KING above", king_threshold,
                 "(het_threshold is used only with called genotypes)\n"))
    }
  }

  # DO THE JOB
  tab <- replicates.report$table.all
  flagged <- .king_flag(tab, king_threshold, het_threshold, loc_threshold,
                        use.het)
  grp <- .king_groups(tab$ind1, tab$ind2, flagged, replicates.report$ind.miss)
  .king_chain_warn(tab, grp, verbose)

  drop_list <- grp$drop
  absent <- setdiff(drop_list, indNames(x))
  if (length(absent) > 0) {
    if (verbose >= 1) {
      cat(warn("  Warning: samples to remove that are not in x (was the report made from another object?):",
               paste(absent, collapse = ", "), "\n"))
    }
    drop_list <- setdiff(drop_list, absent)
  }

  if (length(drop_list) == 0) {
    if (verbose >= 1) {
      cat(report("  No duplicate pairs pass the thresholds; no individuals removed\n"))
      cat(report("Completed:", funname, "\n"))
    }
    return(invisible(x))
  }

  if (verbose >= 2) {
    cat(report("  Removing duplicate samples:", paste(drop_list, collapse = ", "),
               "\n"))
  }

  hold.history <- x@other$history
  xx <- gl.drop.ind(x,
                    ind.list = drop_list,
                    recalc = recalc,
                    mono.rm = mono.rm,
                    verbose = 0)

  # REPORT THE RESULTS
  if (verbose >= 3) {
    cat("  Pairs flagged as duplicates:", sum(flagged), "\n")
    cat("  Initial number of individuals:", nInd(x), "\n")
    cat("  Number of individuals removed:", length(drop_list), "\n")
    cat(report("  Number of individuals retained:", nInd(xx), "\n"))
  }

  # ADD TO HISTORY (a single entry for this call)
  xx@other$history <- hold.history
  nh <- length(xx@other$history)
  xx@other$history[[nh + 1]] <- match.call()

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(invisible(xx))
}
