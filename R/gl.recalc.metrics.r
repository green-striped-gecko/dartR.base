#' @name gl.recalc.metrics
#' @title Recalculates locus metrics when individuals or populations are deleted from a
#'  genlight \{adegenet\} object
#' @family environment

#' @description
#' When individuals, or populations, are deleted from a genlight object, the
#' locus metrics no longer apply. For example, the Call Rate may be different
#' considering the subset of individuals, compared with the full set. This
#' script recalculates those affected locus metrics. For SNP data these are
#' AvgPIC, CallRate, FreqHets, FreqHomRef, FreqHomSnp, OneRatioRef,
#' OneRatioSnp, PICRef, PICSnp and maf. For SilicoDArT (tag presence/absence)
#' data they are CallRate, OneRatio and PIC.

#' @details
#' The script optionally removes resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs).

#' The script returns a genlight object with the recalculated locus metadata.

#' Metrics that are not recalculated fall into two groups. Some are unaffected
#' by the removal of individuals and are correctly left as they are: RepAvg,
#' TrimmedSequence, AlleleSequence, SNP, SnpPosition, clone and uid for SNP
#' data; Qpmr and Reproducibility for SilicoDArT data. Others are affected by
#' the removal of individuals but cannot be recovered from the genotypes,
#' because they derive from raw read counts that a genlight object does not
#' carry: rdepth, AvgCountRef and AvgCountSnp for SNP data; AvgReadDepth and
#' StDevReadDepth for SilicoDArT data. These keep the values read from the DArT
#' file, and can therefore be inconsistent with the recalculated metrics.

#' Precondition on the locus metrics: the function requires
#' \code{x@other$loc.metrics} to be a data frame with one row per locus, the
#' structure that gl.read.dart and gl.compliance.check produce. If the slot is
#' absent or is not a data frame, as in a genlight object assembled by hand, an
#' empty data frame with one row per locus is created to receive the metrics
#' and a message states that the DArT metadata is not present. If the slot is a
#' data frame whose row count does not match the number of loci, the object is
#' not repairable here and the function stops.

#' Monomorphic loci: the monomorphs flag records whether the object is known to
#' be free of monomorphic loci. With mono.rm = FALSE the flag is set from a
#' check made on the metrics just recalculated, so it never reports a state
#' that was not examined. With mono.rm = TRUE the loci are removed and the flag
#' is set by gl.filter.monomorphs; if every locus is monomorphic or scored all
#' NA, none can be removed, and the function reports that condition and
#' completes with the flag FALSE.

#' Messages: the utils.recalc.* helpers are called with verbose = 0 and this
#' function reports on their behalf, so one call produces one set of progress
#' messages rather than one set per helper.

#' History: gl.recalc.metrics is an implementation step of many other
#' dartRverse functions. An entry is appended to \code{x@other$history} only
#' when the function is called directly. A call made from inside another
#' dartRverse function appends nothing, so the history records the calls the
#' user made rather than the steps beneath them.

#' @param x Name of the genlight object containing SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param mono.rm If TRUE, removes monomorphic loci [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log; 3, progress and results summary; 5, full report
#'  [default NULL, adopting the global verbosity set by gl.set.verbosity(), or
#'  2 if no global is set].

#' @return A genlight object with the recalculated locus metadata.

#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#'   gl <- gl.recalc.metrics(testset.gl, verbose=2)

#' @seealso \code{\link{gl.filter.monomorphs}}

#' @export

gl.recalc.metrics <- function(x,
                              mono.rm = FALSE,
                              verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     verbose = verbose)

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)

    # FUNCTION SPECIFIC ERROR CHECKING

    if (!is.logical(mono.rm) ||
        length(mono.rm) != 1 || is.na(mono.rm)) {
        stop(error(
            paste0("Fatal Error: mono.rm must be a single logical ",
                   "value, TRUE or FALSE\n")
        ))
    }

    # The utils.recalc.* helpers read the locus metrics table with $, which
    # partial-matches to loc.metrics.flags when loc.metrics is absent. That
    # writes the metrics into the flags table and returns it as the locus
    # metrics, losing every DArT metadata column without a message. Index the
    # slot exactly here and guarantee a conforming table before the helpers
    # see the object.
    lm <- x@other[["loc.metrics"]]
    if (is.null(lm) || !is.data.frame(lm)) {
        if (verbose >= 1) {
            cat(
                warn(
                    "  Warning: no locus metrics data frame found;",
                    "creating one with",
                    nLoc(x),
                    "rows\n"
                )
            )
            cat(
                warn(
                    "    DArT metadata (AlleleID, TrimmedSequence, rdepth,",
                    "AvgCountRef, AvgCountSnp, RepAvg) is not present in this",
                    "object and cannot be reconstructed from the genotypes\n"
                )
            )
        }
        x@other$loc.metrics <- data.frame(row.names = seq_len(nLoc(x)))
    } else if (nrow(lm) != nLoc(x)) {
        stop(error(
            paste0(
                "Fatal Error: the locus metrics data frame has ",
                nrow(lm),
                " rows for ",
                nLoc(x),
                " loci. Locus metrics must track loci one for one. Run ",
                "gl.compliance.check() on the object first.\n"
            )
        ))
    }

    # The helpers read x@other$loc.metrics.flags$monomorphs before testing it,
    # so an object with no flags list fails on the same path. Seed the flag as
    # FALSE, the "not checked" state; the monomorph check below sets its true
    # value.
    if (is.null(x@other[["loc.metrics.flags"]][["monomorphs"]])) {
        x@other$loc.metrics.flags$monomorphs <- FALSE
    }

    # DO THE JOB

    # Recalculate statistics. The helpers are called with verbose = 0 and this
    # function reports on their behalf, so one call produces one set of
    # messages rather than one set per helper.

    if (datatype == "SNP") {
        x <- utils.recalc.avgpic(x, verbose = 0)
        x <- utils.recalc.callrate(x, verbose = 0)
        x <- utils.recalc.maf(x, verbose = 0)
    }
    if (datatype == "SilicoDArT") {
        x <- utils.recalc.avgpic(x, verbose = 0)
        x <- utils.recalc.callrate(x, verbose = 0)
    }

    if (verbose >= 2) {
        cat(report("  Locus metrics recalculated\n"))
    }

    # Read the refreshed metrics back once, for the reports below and for the
    # monomorph check.
    lm <- x@other[["loc.metrics"]]

    n.allna <- sum(lm$CallRate == 0, na.rm = TRUE)
    if (n.allna > 0 && verbose >= 2) {
        cat(warn(
            "  Warning:",
            n.allna,
            "loci have no scores at all; their metrics are NaN\n"
        ))
    }

    # Monomorphic loci, by the definition gl.filter.monomorphs applies: every
    # scored genotype the same, or the locus scored all NA. Read from the
    # metrics just recalculated, so the genotypes are not scanned again.
    if (datatype == "SNP") {
        mono <- (lm$CallRate == 0) |
            (!is.na(lm$FreqHomRef) & lm$FreqHomRef == 1) |
            (!is.na(lm$FreqHomSnp) & lm$FreqHomSnp == 1)
    } else {
        mono <- (lm$CallRate == 0) |
            (!is.na(lm$OneRatio) &
                 (lm$OneRatio == 0 | lm$OneRatio == 1))
    }
    n.mono <- sum(mono, na.rm = TRUE)

    if (mono.rm) {
        hold.history <- x@other$history
        # gl.filter.monomorphs cannot return a zero-locus object: with every
        # locus monomorphic the subset in gl.drop.loc errors ("Subsetting
        # resulted in zero loci"). Catch that and report the condition.
        x2 <- tryCatch(
            gl.filter.monomorphs(x, verbose = 0),
            error = function(e)
                NULL
        )
        if (is.null(x2)) {
            if (verbose >= 1) {
                cat(warn(
                    "  Warning: all",
                    nLoc(x),
                    "loci are monomorphic or scored all NA; none removed\n"
                ))
            }
            x@other$loc.metrics.flags$monomorphs <- FALSE
        } else {
            x <- x2
            # One entry per user call: discard the entry gl.filter.monomorphs
            # appends for an internal step the user did not invoke.
            x@other$history <- hold.history
            if (verbose >= 2) {
                cat(report("  Monomorphic loci deleted\n"))
            }
        }
    } else {
        # Set the flag from the check just made rather than passing through a
        # value that predates the recalculation.
        x@other$loc.metrics.flags$monomorphs <- (n.mono == 0)
        if (verbose >= 2) {
            if (n.mono > 0) {
                cat(
                    warn(
                        "  Warning:",
                        n.mono,
                        "monomorphic loci (or loci scored all NA) are present",
                        "and retained; use mono.rm = TRUE to remove them\n"
                    )
                )
            } else {
                cat(report("  No monomorphic loci detected\n"))
            }
        }
    }

    # ADD TO HISTORY
    # Only a direct call is recorded. gl.recalc.metrics is an implementation
    # step of many other dartRverse functions; an entry appended on their
    # behalf records a call the user never made, and nested calls multiply
    # entries. Suppress the append when any frame above this one belongs to a
    # dartRverse namespace.
    internal.call <- FALSE
    nf <- sys.nframe()
    if (nf > 1) {
        for (i in seq_len(nf - 1)) {
            env <- tryCatch(environment(sys.function(i)),
                            error = function(e)
                                NULL)
            if (!is.null(env) && isNamespace(env) &&
                grepl("^dartR", environmentName(env))) {
                internal.call <- TRUE
                break
            }
        }
    }
    if (!internal.call) {
        nh <- length(x@other$history)
        x@other$history[[nh + 1]] <- match.call()
    }

    # FLAG SCRIPT END

    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }

    return(x)

}
