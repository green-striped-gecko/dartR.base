#' @name gl.filter.heterozygosity
#' @title Filters individuals with average heterozygosity greater than a
#'  specified upper threshold or less than a specified lower threshold
#' @family matched filter

#' @description
#' Calculates the observed heterozygosity for each individual in a genlight
#' object and filters individuals based on specified threshold values.
#' Use gl.report.heterozygosity to determine the appropriate thresholds.

#' @details
#' Individuals with observed heterozygosity in the interval
#' [t.lower, t.upper] (boundaries included) are retained. Individuals whose
#' heterozygosity cannot be computed (all genotypes missing) are removed,
#' because they cannot be assessed against the thresholds.
#'
#' Removing individuals invalidates the locus metrics that depend on the
#' composition of individuals (CallRate, allele frequencies, PIC values and
#' the like); when individuals are removed, the corresponding locus-metric
#' flags are reset so that downstream functions know to recalculate them
#' (use gl.recalc.metrics to recalculate immediately).

#' @param x A genlight object containing the SNP genotypes [required].
#' @param t.upper Filter individuals > the threshold [default 0.7].
#' @param t.lower Filter individuals < the threshold [default 0].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @return The filtered genlight object.

#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#'  result <- gl.filter.heterozygosity(testset.gl,t.upper=0.06,verbose=3)
#'  tmp <- gl.report.heterozygosity(result,method='ind')

#' @seealso \code{\link{gl.report.heterozygosity}}

#' @importFrom plyr join
#' @export

gl.filter.heterozygosity <- function(x,
                                     t.upper = 0.7,
                                     t.lower = 0,
                                     verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)

    # CHECK DATATYPE
    datatype <-
        utils.check.datatype(x, accept = "SNP", verbose = verbose)

    # SCRIPT SPECIFIC ERROR CHECKING

    if (t.upper < 0 | t.upper > 1) {
        stop(error(
            "Fatal Error: Parameter 't.upper' must lie between 0 and 1\n"
        ))
    }

    if (t.lower < 0 | t.lower > 1) {
        stop(error(
            "Fatal Error: Parameter 't.lower' must lie between 0 and 1\n"
        ))
    }

    if (t.upper < t.lower) {
        if (verbose >= 1) {
            cat(warn(
                "  Warning: Parameter 't.upper' must be greater than parameter
                't.lower', swapping\n"
            ))
        }
        tmp <- t.upper
        t.upper <- t.lower
        t.lower <- tmp
    }

    # Check for monomorphic loci

    if (isFALSE(x@other$loc.metrics.flags$monomorphs)) {
        if (verbose >= 2) {
            cat(
                warn(
                    "  Warning: genlight object contains monomorphic loci which will
                be factored into heterozygosity estimates\n"
                )
            )
        }
    }

    # DO THE JOB

    # Convert to matrix
    m <- as.matrix(x)

    # For each individual determine counts of hets, homs and NAs
    c.na <- array(NA, nInd(x))
    c.hets <- array(NA, nInd(x))
    c.hom0 <- array(NA, nInd(x))
    c.hom2 <- array(NA, nInd(x))
    for (i in 1:nInd(x)) {
        c.na[i] <- sum(is.na(m[i,]))
        c.hets[i] <- sum(m[i,] == 1, na.rm = TRUE) / (nLoc(x) - c.na[i])
        c.hom0[i] <- sum(m[i,] == 0, na.rm = TRUE) / (nLoc(x) - c.na[i])
        c.hom2[i] <- sum(m[i,] == 2, na.rm = TRUE) / (nLoc(x) - c.na[i])
    }

    # Individuals with all genotypes missing have undefined heterozygosity
    # (NaN) and cannot be assessed against the thresholds; they are removed
    not.assessable <- !is.finite(c.hets)

    if (verbose >= 2) {
        cat(
            report(
                "  Retaining individuals with heterozygosity in the range",
                t.lower,
                "to",
                t.upper,
                "\n"
            )
        )
        cat(paste(
            "  Minimum individual heterozygosity",
            round(min(c.hets, na.rm = TRUE), 4),
            "\n"
        ))
        cat(paste(
            "  Maximum individual heterozygosity",
            round(max(c.hets, na.rm = TRUE), 4),
            "\n"
        ))

    }
    keep <- !not.assessable & c.hets >= t.lower & c.hets <= t.upper
    x.kept <- x[keep,]
    upper.out <- !not.assessable & c.hets > t.upper
    lower.out <- !not.assessable & c.hets < t.lower
    if (any(upper.out)) {
        x.discarded.upper <- x[upper.out,]
    }
    if (any(lower.out)) {
        x.discarded.lower <- x[lower.out,]
    }

    # REPORT THE RESULTS
    if (verbose >= 3) {
        cat("  Initial number of individuals:", nInd(x), "\n")

        if (any(upper.out)) {
            cat(
                "  Number of outlier individuals (heterozygosity  >",
                t.upper,
                "):",
                nInd(x.discarded.upper),
                "\n"
            )
            cat("    Deleted:",
                paste(paste0(
                    indNames(x.discarded.upper),
                    "[",
                    as.character(pop(x.discarded.upper)),
                    "]"
                ), collapse = ", "),
                "\n")
        } else {
            if (!(t.upper == 1)) {
                cat("  Zero outlier individuals with heterozygosity >",
                    t.upper,
                    "\n")
            }
        }

        if (any(lower.out)) {
            cat(
                "  Number of outlier individuals:",
                nInd(x.discarded.lower),
                "with heterozygosity <",
                t.lower,
                "\n"
            )
            cat("    Deleted:",
                paste(paste0(
                    indNames(x.discarded.lower),
                    "[",
                    as.character(pop(x.discarded.lower)),
                    "]"
                ), collapse = ", "),
                "\n")
        } else {
            if ((!t.lower == 0)) {
                cat("  No outlying individuals with heterozygosity <",
                    t.lower,
                    "\n")
            }
        }
        if (any(not.assessable)) {
            cat("  Number of individuals removed for undefined heterozygosity
          (all genotypes missing):",
                sum(not.assessable),
                "\n")
        }
        cat(report("  Number of individuals retained:", nInd(x.kept), "\n"))
    }

    # RESET THE FLAGS for metrics that depend on the composition of
    # individuals, when individuals have been removed
    if (nInd(x.kept) != nInd(x)) {
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
    }

    # ADD ACTION TO HISTORY

    nh <- length(x.kept@other$history)
    x.kept@other$history[[nh + 1]] <- match.call()

    # FLAG SCRIPT END

    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }

    return(invisible(x.kept))
}
