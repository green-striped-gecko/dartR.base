#' @name gl.sample
#' @title Samples individuals from populations
#' @family data manipulation

#' @description
#' Draws a subsample of individuals from a genlight object, with or without
#' replacement, either from each population separately or from the object as a
#' whole.

#' @details
#' This function supports a bootstrap approach in dartR. For a bootstrap it is
#' often desirable to draw a defined number of individuals from each of the
#' populations in a genlight object, calculate a quantity for that subset, and
#' repeat the draw many times.
#'
#' Parameter nsample has two meanings, set by onepop:
#' \itemize{
#'  \item onepop = FALSE: nsample is a per-population count. The result holds
#'  nsample individuals from every population, including a population with only
#'  one member, so nInd is nsample times the number of populations.
#'  \item onepop = TRUE: population structure is ignored and nsample is a
#'  whole-object total. The result holds nsample individuals drawn from all
#'  individuals, each carrying its original population label.
#' }
#' The default nsample follows the same split: the size of the smallest
#' population when onepop = FALSE, and nInd(x) when onepop = TRUE. Population
#' factor levels with no members are ignored.
#'
#' An individual drawn more than once under replace = TRUE would give the result
#' duplicate individual names, and duplicate names corrupt any downstream join
#' on individual identity. Every returned individual is therefore renamed with a
#' zero-padded ordinal prefix recording its position in the sample, for example
#' 01_AA013220, which makes the names unique by construction. The same names are
#' written to @@other$ind.metrics$id so that the two agree.
#'
#' Resampling individuals invalidates every locus metric computed across
#' individuals, so all entries of @@other$loc.metrics.flags are set to FALSE.
#' The stored locus metrics themselves are left as they are, which keeps the
#' cost of a bootstrap replicate down; downstream report and filter functions
#' recalculate what they need because the flags tell them the stored values are
#' stale. Call gl.recalc.metrics() to refresh the metrics in place.

#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param nsample Number of individuals to draw, per population when
#' onepop = FALSE and in total when onepop = TRUE. Must be a positive whole
#' number [default NULL, the size of the smallest population when
#' onepop = FALSE, nInd(x) when onepop = TRUE].
#' @param replace If TRUE, sampling is with replacement [default TRUE].
#' @param onepop If TRUE, ignore the population assignments of the genlight
#' object and draw from all individuals as a single pool [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @return A genlight object holding the drawn individuals, in draw order.
#' Individual-indexed metadata in @@other (ind.metrics, latlon and the like) is
#' subset to the drawn individuals, repeating rows for individuals drawn more
#' than once; @@other$loc.metrics is carried over unchanged because the loci are
#' unchanged; all entries of @@other$loc.metrics.flags are set to FALSE; and the
#' call is appended to @@other$history. Individuals are renamed with a
#' zero-padded ordinal prefix (see Details).

#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' \donttest{
#' # bootstrap for 2 possums populations to check effect of sample size on
#' # fixed alleles
#' gl.set.verbosity(0)
#' if (isTRUE(getOption("dartR_fbm"))) possums.gl <- gl.gen2fbm(possums.gl)
#' pp <- possums.gl[c(1:30,91:120),]
#' nrep <- 1:10
#' nss <- seq(1,10,2)
#' res <- expand.grid(nrep=nrep, nss=nss)
#' for (i in 1:nrow(res)) {
#' dummy <- gl.sample(pp, nsample=res$nss[i], replace=TRUE)
#' pas <- gl.report.pa(dummy, plot.display= FALSE)
#' res$fixed[i] <- pas$fixed[1]
#' }
#' boxplot(fixed ~ nss, data=res)
#'}
#' @export

gl.sample <- function(x,
                      nsample = NULL,
                      replace = TRUE,
                      onepop = FALSE,
                      verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname, verbose = verbose)

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)

    # FUNCTION SPECIFIC ERROR CHECKING

    if (length(onepop) != 1 || !is.logical(onepop) || is.na(onepop)) {
        stop(error("Fatal Error: onepop must be TRUE or FALSE\n"))
    }
    if (length(replace) != 1 || !is.logical(replace) || is.na(replace)) {
        stop(error("Fatal Error: replace must be TRUE or FALSE\n"))
    }

    # Population sizes, counting only levels that have members. An unused factor
    # level is not a population, so it neither sets the default nsample nor
    # contributes a draw.
    pop.sizes <- table(pop(x))
    pop.sizes <- pop.sizes[pop.sizes > 0]
    if (!onepop && length(pop.sizes) == 0) {
        stop(error(
            "Fatal Error: the genlight object has no population assignments;",
            "use onepop = TRUE or assign populations first\n"
        ))
    }

    # Resolve the default deliberately, before nsample is used. The former
    # signature default min(table(pop(x))) was a lazily evaluated promise that
    # resolved only after pop(x) had been overwritten under onepop = TRUE, so
    # the default depended on evaluation order rather than on intent.
    if (is.null(nsample)) {
        nsample <- if (onepop) nInd(x) else min(pop.sizes)
    }

    if (length(nsample) != 1 || !is.numeric(nsample) || is.na(nsample)) {
        stop(error(
            "Fatal Error: nsample must be a single positive whole number\n"
        ))
    }
    if (nsample != round(nsample)) {
        stop(error(
            "Fatal Error: nsample must be a whole number, it is a count of",
            "individuals and not a proportion. Value provided:",
            nsample,
            "\n"
        ))
    }
    nsample <- as.integer(round(nsample))
    if (nsample < 1) {
        stop(error(
            "Fatal Error: nsample must be 1 or more. Value provided:",
            nsample,
            "\n"
        ))
    }
    if (!replace) {
        if (onepop) {
            if (nsample > nInd(x)) {
                stop(error(
                    "Fatal Error: nsample",
                    paste0("(", nsample, ")"),
                    "exceeds the number of individuals",
                    paste0("(", nInd(x), ")"),
                    "and replace is FALSE. Reduce nsample or set",
                    "replace = TRUE\n"
                ))
            }
        } else if (nsample > min(pop.sizes)) {
            smallest <- names(pop.sizes)[which.min(pop.sizes)]
            stop(error(
                "Fatal Error: nsample",
                paste0("(", nsample, ")"),
                "exceeds the size of the smallest population,",
                smallest,
                "with",
                min(pop.sizes),
                "individuals, and replace is FALSE. Reduce nsample or set",
                "replace = TRUE\n"
            ))
        }
    }

    # DO THE JOB

    # Index vectors of the individuals available to each draw: one pool per
    # population, or a single pool covering the whole object under onepop.
    if (onepop) {
        pools <- list(seq_len(nInd(x)))
    } else {
        pops <- pop(x)
        pools <- lapply(names(pop.sizes), function(p) which(pops == p))
    }

    # Draw from each pool. sample.int(length(idx), ...) indexes into the pool
    # rather than handing the pool to sample(), because sample(v, n) treats a
    # length-1 v as the range 1:v; that is what replaced a single-member
    # population with unrelated individuals drawn from elsewhere in the object.
    samps <- unlist(lapply(pools, function(idx) {
        idx[sample.int(length(idx), nsample, replace = replace)]
    }))

    # Subset with a positive index vector, in draw order. The dartR "[" method
    # carries @other with it, so individual-indexed metadata (ind.metrics,
    # latlon, and anything else with one row or element per individual) tracks
    # the drawn individuals, repeated rows included.
    #
    # "[" cannot return more individuals than the object holds: it assigns
    # ploidy through the accessor after @ind.names has already been replaced, so
    # an index vector longer than nInd(x) fails with "'names' attribute [n] must
    # be the same length as the vector [nInd]". A draw that stays within nInd(x)
    # -- every draw of nsample <= the smallest population, and so every ordinary
    # bootstrap -- therefore goes through "[" untouched. A larger draw is taken
    # as the distinct individuals drawn, then expanded to the full draw by
    # repeating rows, which gives the same result by the same indices.
    if (length(samps) <= nInd(x)) {
        xx <- x[samps, ]
    } else {
        drawn <- unique(samps)
        pos <- match(samps, drawn)
        xx <- x[drawn, ]
        if (.has_fbm(xx)) {
            xx@fbm <- bigstatsr::big_copy(
                xx@fbm,
                ind.row = pos,
                ind.col = seq_len(nLoc(xx)),
                backingfile = tempfile("geno_")
            )
        } else {
            xx@gen <- xx@gen[pos]
        }
        xx@ind.names <- xx@ind.names[pos]
        if (!is.null(xx@ploidy)) {
            xx@ploidy <- xx@ploidy[pos]
        }
        if (!is.null(xx@pop)) {
            xx@pop <- factor(xx@pop[pos], levels = levels(xx@pop))
        }
        if (!is.null(xx@strata)) {
            xx@strata <- xx@strata[pos, , drop = FALSE]
        }
        # Expand the individual-indexed elements of @other to match. The slots
        # that are not indexed by individual are named rather than inferred from
        # their shape, so a coincidence of dimensions cannot scramble them.
        n.drawn <- length(drawn)
        not.ind <- c("loc.metrics", "loc.metrics.flags", "history", "verbose")
        expand <- !(names(xx@other) %in% not.ind)
        xx@other[expand] <- lapply(xx@other[expand], function(obj) {
            if (!is.null(dim(obj)) && nrow(obj) == n.drawn) {
                obj[pos, , drop = FALSE]
            } else if (is.null(dim(obj)) && length(obj) == n.drawn) {
                obj <- obj[pos]
                if (is.factor(obj)) factor(obj) else obj
            } else {
                obj
            }
        })
    }

    # The loci are untouched, so restore loc.metrics from the input verbatim
    # rather than relying on the subsetting method to have left it alone.
    if (!is.null(x@other$loc.metrics)) {
        xx@other$loc.metrics <- x@other$loc.metrics
    }

    # The individuals have changed, so every locus metric computed across
    # individuals is now stale. Flag them all rather than pay for recalculation
    # on every bootstrap replicate.
    if (!is.null(xx@other$loc.metrics.flags)) {
        xx@other$loc.metrics.flags[] <-
            lapply(xx@other$loc.metrics.flags, function(f) FALSE)
    }

    # Rename the drawn individuals with a zero-padded ordinal prefix. This is
    # what keeps the names unique when an individual is drawn more than once,
    # and ind.metrics$id is written through so that identity joins still work.
    n10 <- nchar(as.character(nInd(xx)))
    lzs <- paste0("%0", as.character(n10), "d")
    new.names <- paste0(sprintf(lzs, seq_len(nInd(xx))), "_", indNames(xx))
    indNames(xx) <- new.names
    if (!is.null(xx@other$ind.metrics) &&
        "id" %in% names(xx@other$ind.metrics)) {
        xx@other$ind.metrics$id <- new.names
    }

    if (verbose >= 2) {
        if (onepop) {
            cat(report(
                "  Drew",
                nsample,
                "individuals from all",
                nInd(x),
                "individuals,",
                if (replace) "with" else "without",
                "replacement, ignoring population structure\n"
            ))
        } else {
            cat(report(
                "  Drew",
                nsample,
                "individuals from each of",
                length(pools),
                "populations,",
                if (replace) "with" else "without",
                "replacement\n"
            ))
        }
        cat(report(
            "  Locus metrics flags set to FALSE; run gl.recalc.metrics() to",
            "refresh the locus metrics\n"
        ))
    }

    # ADD TO HISTORY
    nh <- length(xx@other$history)
    xx@other$history[[nh + 1]] <- match.call()

    # FLAG SCRIPT END
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }

    return(xx)
}
