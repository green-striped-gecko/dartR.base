#' @name gl.filter.contamination
#' @title Removes individuals flagged as cross-contaminated samples
#' @family matched filter
#'
#' @description
#' Runs the contamination screen of \code{\link{gl.report.contamination}}
#' and removes the individuals it flags. Use the report first to check the
#' candidate list and tune the thresholds; the filter takes the same
#' parameters.
#'
#' @details
#' The screen is described in \code{\link{gl.report.contamination}}. In
#' brief, an individual is a "suspect" when its heterozygosity is a robust
#' outlier within its population and exceeds the population median by
#' more than \code{min.excess}; a suspect whose strongest residual-kinship
#' partner sits in an adjacent plate well is "adjacent"; an individual whose
#' rare-allele burden is elevated but whose heterozygosity is normal is
#' "rare-only", which points to population structure or mislabelling rather
#' than contamination.
#'
#' \code{flag} names the classes to remove. The default removes "suspect"
#' and "adjacent", the two classes with the heterozygosity signal that
#' contamination always produces. Add "rare-only" to also remove individuals
#' that carry an excess of alleles rare in their assigned population, or set
#' \code{flag = "adjacent"} to remove only the suspects with plate evidence.
#'
#' The tier 3 depth \code{pattern} (dose, flat, saturated, unclear) is
#' printed for each removed individual but does not change the flag or the
#' filter: a "flat" pattern points to a hybrid, admixed or mislabelled
#' animal rather than a contaminant, so review such cases before removing
#' them, for example by dropping them from the screen with
#' \code{\link{gl.keep.ind}} or adding them back afterwards.
#'
#' The screen is a candidate list, not a verdict: a population whose
#' members differ widely in heterozygosity (inbreeding, admixture,
#' relatives) can yield suspects that are not contaminated. Check the
#' report before filtering.
#'
#' Removing individuals invalidates the locus metrics that depend on the
#' composition of individuals (CallRate, allele frequencies, PIC values and
#' the like). When individuals are removed the corresponding locus-metric
#' flags are reset so that downstream functions know to recalculate them,
#' or the metrics are recalculated immediately with \code{recalc = TRUE}.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param flag Classes of individual to remove, any of "suspect",
#' "adjacent" and "rare-only" [default c("suspect", "adjacent")].
#' @param rare.freq Allele frequency in the rest of the population below
#' which an allele counts as rare; must lie in (0, 0.5) [default 0.02].
#' @param min.n Minimum number of other individuals in the population called
#' at a locus for that locus to enter the rare-allele burden; at least 2
#' [default 5].
#' @param z.flag Robust z-score above which a statistic is an outlier
#' [default 3].
#' @param min.excess Absolute amount by which heterozygosity must also exceed
#' the population median to be flagged [default 0.02].
#' @param rare.min.excess The same absolute margin for the rare-allele
#' burden [default 0.005].
#' @param min.share Minimum number of rare-allele loci an individual needs
#' before its source is chosen by allele sharing rather than by residual
#' kinship alone [default 20].
#' @param share.tol Sharing scores within this distance of the best are
#' treated as ties and separated by residual kinship [default 0.05].
#' @param depth.ratio Foreign-allele rate in the deepest quartile of loci
#' over the shallowest, at or above which the depth pattern is "dose"
#' [default 1.5].
#' @param depth.flat The same ratio below which the pattern is "flat"
#' [default 1.25].
#' @param plate Data frame with columns id, plate and well giving the plate
#' position of each individual; when NULL, positions are taken from
#' \code{@@other$ind.metrics} if present [default NULL].
#' @param recalc Recalculate the locus metadata statistics if any individuals
#' are deleted [default FALSE].
#' @param mono.rm Remove monomorphic loci after individuals are deleted
#' [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @return The genlight object with the flagged individuals removed.
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' out <- gl.report.contamination(testset.gl, plot.display = FALSE)
#' gl <- gl.filter.contamination(testset.gl, verbose = 3)
#' nInd(testset.gl) - nInd(gl)
#'
#' @seealso \code{\link{gl.report.contamination}},
#' \code{\link{gl.filter.heterozygosity}}
#' @export

gl.filter.contamination <- function(x,
                                    flag = c("suspect", "adjacent"),
                                    rare.freq = 0.02,
                                    min.n = 5,
                                    z.flag = 3,
                                    min.excess = 0.02,
                                    rare.min.excess = 0.005,
                                    min.share = 20,
                                    share.tol = 0.05,
                                    depth.ratio = 1.5,
                                    depth.flat = 1.25,
                                    plate = NULL,
                                    recalc = FALSE,
                                    mono.rm = FALSE,
                                    verbose = NULL) {

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  classes <- c("suspect", "adjacent", "rare-only")
  if (length(flag) == 0 || !all(flag %in% classes)) {
    stop(error("  flag must be one or more of ",
               paste(shQuote(classes), collapse = ", "), "\n"))
  }

  # DO THE JOB
  # The screen validates the remaining parameters and reports at verbose 0
  # so that this function controls what is printed
  screen <- gl.report.contamination(x,
                                    rare.freq = rare.freq,
                                    min.n = min.n,
                                    z.flag = z.flag,
                                    min.excess = min.excess,
                                    rare.min.excess = rare.min.excess,
                                    min.share = min.share,
                                    share.tol = share.tol,
                                    depth.ratio = depth.ratio,
                                    depth.flat = depth.flat,
                                    plate = plate,
                                    plot.display = FALSE,
                                    verbose = 0)
  ind <- screen$ind
  drop <- ind[ind$flag %in% flag, , drop = FALSE]

  # The screen ran silently, so repeat the caveats that limit what can be
  # removed
  small <- names(which(table(pop(x)) < 2))
  if (length(small) > 0 && verbose >= 2) {
    cat(warn("  Warning: populations with a single individual are not",
             "screened and none of their members is removed:",
             paste(small, collapse = ", "), "\n"))
  }
  # Without plate positions no individual can be "adjacent"; this empties
  # the result when "adjacent" is the only class requested
  if ("adjacent" %in% flag && all(is.na(ind$adjacent))) {
    lvl <- if ("suspect" %in% flag) 2 else 1
    if (verbose >= lvl) {
      cat(warn("  Warning: no plate positions found, so no individual can",
               "be flagged \"adjacent\"\n"))
    }
  }
  keep <- !(indNames(x) %in% drop$id)

  if (verbose >= 2) {
    cat(report("  Removing individuals flagged as",
               paste(flag, collapse = ", "), "\n"))
  }

  x.kept <- x[keep, ]

  if (any(!keep)) {
    if (mono.rm) {
      x.kept <- gl.filter.monomorphs(x.kept, verbose = 0)
    }
    if (recalc) {
      x.kept <- gl.recalc.metrics(x.kept, verbose = 0)
    } else {
      # Reset the flags for metrics that depend on the composition of
      # individuals
      x.kept@other$loc.metrics.flags$AvgPIC <- FALSE
      x.kept@other$loc.metrics.flags$OneRatioRef <- FALSE
      x.kept@other$loc.metrics.flags$OneRatioSnp <- FALSE
      x.kept@other$loc.metrics.flags$PICRef <- FALSE
      x.kept@other$loc.metrics.flags$PICSnp <- FALSE
      x.kept@other$loc.metrics.flags$CallRate <- FALSE
      x.kept@other$loc.metrics.flags$maf <- FALSE
      x.kept@other$loc.metrics.flags$FreqHets <- FALSE
      x.kept@other$loc.metrics.flags$FreqHomRef <- FALSE
      x.kept@other$loc.metrics.flags$FreqHomSnp <- FALSE
      x.kept@other$loc.metrics.flags$allna <- FALSE
    }
  }

  # REPORT THE RESULTS
  if (verbose >= 3) {
    cat("  Initial number of individuals:", nInd(x), "\n")
    cat("  Number of individuals removed:", nrow(drop), "\n")
    if (nrow(drop) > 0) {
      show <- c("id", "pop", "well", "het", "het.excess", "rare.burden",
                "top.partner", "partner.well", "pattern", "flag")
      print(drop[, show], row.names = FALSE)
    }
    if (mono.rm) {
      cat("  Number of loci retained after removing monomorphs:",
          nLoc(x.kept), "\n")
    }
    cat(report("  Number of individuals retained:", nInd(x.kept), "\n"))
  }

  # ADD ACTION TO HISTORY
  nh <- length(x.kept@other$history)
  x.kept@other$history[[nh + 1]] <- match.call()

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(invisible(x.kept))
}
