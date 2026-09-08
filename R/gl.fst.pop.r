#' @name gl.fst.pop
#' @title Calculates pairwise Fst between populations in a genlight object
#' @family distance

#' @description
#' Calculates the pairwise fixation index Fst between every pair of
#' populations in a genlight object, with optional bootstrap confidence
#' intervals and p-values obtained by resampling loci.

#' @param x Name of the genlight containing the SNP genotypes [required].
#' @param nboots Number of bootstrap replicates over loci used to generate
#' confidence intervals and p-values. A whole number; 1 returns the point
#' estimates alone [default 1].
#' @param percent Percentile to calculate the confidence interval around,
#' greater than 0 and less than 100 [default 95].
#' @param nclusters Number of processor threads or cores to use during
#' calculation of the point estimates [default 1].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @details
#' The statistic is Weir and Cockerham's (1984) theta, estimated by
#' \code{\link[StAMPP]{stamppFst}} from allele frequencies, observed
#' heterozygosities and per-locus sample sizes computed separately for each
#' population. Theta is a variance-components estimator that corrects for
#' unequal sample sizes, so it is not the same quantity as the Nei (1987)
#' Gst-family statistics reported by \code{\link{gl.report.fstat}} under the
#' same word "Fst"; the two can differ substantially on unevenly sampled
#' data. Loci missing from one member of a pair contribute nothing to that
#' pair.
#'
#' \strong{Bootstrap.} When \code{nboots} is greater than 1 the function
#' resamples \emph{loci} with replacement, nLoc at a time, and recomputes
#' theta on each replicate with the same estimator that produced the point
#' estimate. The locus is the unit of resampling because it is the unit the
#' variance components are summed over. The replicates are drawn in the
#' calling session, so \code{set.seed()} immediately before the call makes
#' the confidence limits and the p-values reproducible. There is no
#' \code{seed} parameter. The point estimates do not depend on the bootstrap
#' and are the same whatever \code{nboots} is set to.
#'
#' \strong{Confidence limits.} The limits are the order statistics of the
#' sorted replicates at positions \code{ceiling(a * nboots)} and
#' \code{floor((1 - a) * nboots)}, where \code{a = ((100 - percent)/100)/2}.
#' At the default \code{percent = 95} the lower position is 1 for every
#' \code{nboots} up to 40, so with 40 or fewer replicates the reported lower
#' limit is the smallest replicate rather than the 2.5th percentile and the
#' interval is the observed range. The function warns at
#' \code{verbose >= 1} whenever the requested percentile cannot be reached.
#'
#' \strong{P-values.} The p-value for a pair is the one-tailed bootstrap
#' fraction \code{mean(replicates <= 0)}: the proportion of replicate thetas
#' at or below zero. Its resolution is \code{1/nboots}, and it can be
#' exactly 0 or exactly 1. A p of 1 means every replicate was at or below
#' zero, which is what undifferentiated populations produce; it is a failure
#' to reject, not evidence of no differentiation. No correction for multiple
#' testing is applied across the \code{nPop*(nPop-1)/2} pairs tested, which
#' is 435 pairs for a 30-population object. Apply one yourself if the
#' pairwise tests are to be read as a family.
#'
#' \strong{Non-finite pairs.} A pair returns NaN when the estimator cannot
#' be formed, most often because a population holds a single individual or
#' because no locus is polymorphic in the pair. Such pairs are named at
#' \code{verbose >= 1}.

#' @return A base matrix of class \code{c("matrix", "array")} when
#' \code{nboots} is 1: nPop x nPop, with the pairwise theta values in the
#' lower triangle and NA on the diagonal and in the upper triangle. It is
#' not an object of class \code{dist}. Rows and columns are named by
#' population in order of first appearance in the data, not in factor-level
#' order. Use \code{as.matrix(as.dist(fsts))} for a symmetric square matrix.
#'
#' When \code{nboots} is greater than 1, a list of three elements:
#' \itemize{
#' \item \code{Fsts} - the same lower-triangular nPop x nPop matrix.
#' \item \code{Pvalues} - an nPop x nPop matrix of p-values, lower triangle
#' populated.
#' \item \code{Bootstraps} - a data frame with one row per pair of
#' populations and columns "Population1", "Population2", the \code{nboots}
#' replicate thetas sorted in increasing order, "Lower bound CI limit",
#' "Upper bound CI limit", "p-value" and "Fst".
#' }
#' @importFrom StAMPP stamppFst
#' @export
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @references
#' Weir, B.S. and Cockerham, C.C. (1984). Estimating F-statistics for the
#' analysis of population structure. Evolution 38, 1358-1370.

#' @seealso \code{\link{gl.report.fstat}}, which reports the Nei (1987)
#' Gst-family statistics Fst, Fstp, Dest and Gst_H rather than the Weir and
#' Cockerham theta returned here.

#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' test <- gl.filter.callrate(platypus.gl,threshold = 1)
#' test <- gl.filter.monomorphs(test)
#' out <- gl.fst.pop(test, nboots=1)

gl.fst.pop <- function(x,
                       nboots = 1,
                       percent = 95,
                       nclusters = 1,
                       verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # FLAG SCRIPT START [approved F12]
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     verbose = verbose)

    # CHECK DATATYPE
    # theta is computed from allele dosages scaled by ploidy; on
    # presence/absence data the observed-heterozygosity term is identically
    # zero and per-population sample size is halved, so the number returned
    # is not an Fst [approved F2]
    datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

    # FUNCTION SPECIFIC ERROR CHECKING [approved F6]
    if (nPop(x) < 2) {
        stop(error(
            "Fatal Error: pairwise Fst requires at least two populations;",
            "the object holds", nPop(x),
            "- assign populations with gl.reassign.pop or gl.define.pop\n"
        ))
    }

    if (!is.numeric(nboots) || length(nboots) != 1 || is.na(nboots) ||
        nboots < 1 || nboots != round(nboots)) {
        stop(error(
            "Fatal Error: nboots must be a single whole number of 1 or",
            "more; 1 returns the point estimates alone\n"
        ))
    }

    if (!is.numeric(percent) || length(percent) != 1 || is.na(percent) ||
        percent <= 0 || percent >= 100) {
        stop(error(
            "Fatal Error: percent must be a single number greater than 0",
            "and less than 100\n"
        ))
    }

    if (!is.numeric(nclusters) || length(nclusters) != 1 ||
        is.na(nclusters) || nclusters < 1 || nclusters != round(nclusters)) {
        stop(error(
            "Fatal Error: nclusters must be a single whole number of 1 or",
            "more\n"
        ))
    }

    # the lower confidence limit is an order statistic of the sorted
    # replicates; when its position rounds down to 1 the limit is the sample
    # minimum and the interval is the observed range, not a percentile
    # interval [approved F6]
    alpha <- ((100 - percent) / 100) / 2
    if (nboots > 1 && ceiling(alpha * nboots) <= 1 && verbose >= 1) {
        cat(warn(
            "  Warning: with nboots =", nboots, "and percent =", percent,
            "the lower confidence limit is the smallest of the", nboots,
            "replicates, not the", paste0(format(alpha * 100), "th"),
            "percentile, and the reported interval is the observed range.",
            "At least", floor(1 / alpha) + 1,
            "replicates are needed to reach that percentile\n"
        ))
    }

    # DO THE JOB

    #!# intermediate fbm fix
    if (!is.null(.fbm_or_null(x))) x <- gl.fbm2gen(x)

    class(x)<- "genlight" #needs to be genlight due to stampp

    # ---- point estimates ------------------------------------------------
    # StAMPP computes the reported theta. It is called with nboots = 1 so
    # that only the estimator runs; the bootstrap below is this function's
    # own [approved F1]
    fstmat <-
        stamppFst(x,
                  nboots = 1,
                  percent = percent,
                  nclusters = nclusters)

    pops <- rownames(fstmat)
    npops <- length(pops)

    # lower triangle in column-major order: pair k compares population
    # index1[k] (row) with index2[k] (column)
    index1 <- index2 <- integer(0)
    for (j in seq_len(npops - 1)) {
        index1 <- c(index1, (j + 1):npops)
        index2 <- c(index2, rep(j, npops - j))
    }
    npairs <- length(index1)

    if (nboots > 1) {
        # ---- bootstrap over loci ----------------------------------------
        # The bootstrap is computed here rather than by StAMPP, which runs
        # it inside foreach %dopar% on a PSOCK cluster of its own making.
        # Those worker RNG streams are never seeded from the calling
        # session, so set.seed() had no effect on the confidence limits or
        # the p-values and no reported result could be reproduced. Drawing
        # the replicate indices in this session restores reproducibility
        # [approved F1].
        #
        # The unit of resampling is the locus. Per-population allele
        # frequency, observed heterozygosity and sample size are computed
        # once per locus from the DECODED genotype matrix, and a replicate
        # is a resample of the columns of those three matrices. Nothing is
        # subsetted out of the genlight with repeated indices, so adegenet's
        # SNPbin "[" method -- which drops NA when an index repeats, as a
        # draw with replacement always does -- is never involved.

        # per-population, per-locus terms, in the form StAMPP builds them
        stampp.terms <- function(gm, ploidy.vec, pop.vec, pop.order) {
            gm <- gm * (1 / ploidy.vec)  # dosage to allele frequency
            np <- length(pop.order)
            p <- oh <- ninds <-
                matrix(NA_real_, nrow = np, ncol = ncol(gm),
                       dimnames = list(pop.order, NULL))
            w <- ploidy.vec / 2
            for (i in seq_len(np)) {
                sel <- pop.vec == pop.order[i]
                g <- gm[sel, , drop = FALSE]
                wi <- w[sel]
                ninds[i, ] <- colSums((!is.na(g)) * wi, na.rm = TRUE)
                p[i, ] <- colMeans(g, na.rm = TRUE)
                het <- (g * (g <= 0.5)) + ((g - 0.5) * (g > 0.5 & g != 1))
                oh[i, ] <- 2 * colSums(het * wi, na.rm = TRUE)
            }
            list(p = p, oh = oh / ninds, ninds = ninds)
        }

        # Weir & Cockerham (1984) theta for every pair, from those terms
        wc.theta <- function(p, oh, ninds, nloc) {
            r <- 2
            n1 <- ninds[index1, ]
            n2 <- ninds[index2, ]
            p1 <- p[index1, ]
            p2 <- p[index2, ]
            oh1 <- oh[index1, ]
            oh2 <- oh[index2, ]
            n.bar <- (n1 + n2) / r
            nc <- (r * n.bar) - (((n1 ^ 2) + (n2 ^ 2)) / (r * n.bar))
            p.bar <- ((n1 * p1) / (r * n.bar)) + ((n2 * p2) / (r * n.bar))
            s.square <- ((n1 * ((p1 - p.bar) ^ 2)) / n.bar) +
                ((n2 * ((p2 - p.bar) ^ 2)) / n.bar)
            h.bar <- ((n1 * oh1) / (r * n.bar)) + ((n2 * oh2) / (r * n.bar))
            a <- (n.bar / nc) * (s.square - (1 / (n.bar - 1)) *
                ((p.bar * (1 - p.bar)) - (((r - 1) / r) * s.square) -
                     ((1 / 4) * h.bar)))
            b <- (n.bar / (n.bar - 1)) * ((p.bar * (1 - p.bar)) -
                (((r - 1) / r) * s.square) -
                (((2 * n.bar - 1) / (4 * n.bar)) * h.bar))
            cw <- (1 / 2) * h.bar
            bad <- which(!is.finite(a) | !is.finite(b) | !is.finite(cw))
            a[bad] <- NA
            b[bad] <- NA
            cw[bad] <- NA
            if (nloc > 1) {
                if (npops > 2) {
                    rowSums(a, na.rm = TRUE) /
                        (rowSums(a, na.rm = TRUE) +
                             rowSums(b, na.rm = TRUE) +
                             rowSums(cw, na.rm = TRUE))
                } else {
                    sum(a, na.rm = TRUE) /
                        (sum(a, na.rm = TRUE) + sum(b, na.rm = TRUE) +
                             sum(cw, na.rm = TRUE))
                }
            } else {
                a / (a + b + cw)
            }
        }

        gm <- as.matrix(x)
        nloc <- ncol(gm)
        tt <- stampp.terms(gm, ploidy(x), as.character(pop(x)), pops)

        if (verbose >= 2) {
            cat(report("  Bootstrapping", nboots, "replicates over", nloc,
                       "loci\n"))
        }

        reps <- matrix(NA_real_, nrow = npairs, ncol = nboots)
        for (d in seq_len(nboots)) {
            boot.index <- sample(seq_len(nloc), nloc, replace = TRUE)
            reps[, d] <- wc.theta(tt$p[, boot.index, drop = FALSE],
                                  tt$oh[, boot.index, drop = FALSE],
                                  tt$ninds[, boot.index, drop = FALSE],
                                  nloc)
        }

        # sorted replicates, order-statistic limits, and the one-tailed
        # bootstrap fraction at or below zero
        sorted <- t(apply(reps, 1, sort, na.last = TRUE))
        lowerper <- sorted[, ceiling(alpha * nboots)]
        upperper <- sorted[, floor((1 - alpha) * nboots)]
        pval <- rowSums(sorted <= 0, na.rm = TRUE) / nboots

        pvalues <- matrix(NA_real_, nrow = npops, ncol = npops,
                          dimnames = list(pops, pops))
        pvalues[cbind(index1, index2)] <- pval

        boots <- cbind.data.frame(
            pops[index2],
            pops[index1],
            sorted,
            lowerper,
            upperper,
            pval,
            fstmat[cbind(index1, index2)],
            stringsAsFactors = FALSE
        )
        colnames(boots) <- c("Population1", "Population2",
                             as.character(seq_len(nboots)),
                             "Lower bound CI limit", "Upper bound CI limit",
                             "p-value", "Fst")
        rownames(boots) <- NULL

        out <- list(Fsts = fstmat, Pvalues = pvalues, Bootstraps = boots)
    } else {
        out <- fstmat
    }

    # pairs the estimator cannot form come back as NaN; an unannounced NaN
    # cell reads downstream as missing data [approved F7]
    bad.pairs <- which(!is.finite(fstmat[cbind(index1, index2)]))
    if (length(bad.pairs) > 0 && verbose >= 1) {
        pair.txt <- paste0(pops[index1[bad.pairs]], " vs ",
                           pops[index2[bad.pairs]])
        cat(warn(
            "  Warning:", length(bad.pairs), "of", npairs,
            "population pairs returned a non-finite Fst, most likely",
            "because a population holds a single individual or no locus is",
            "polymorphic in the pair:",
            paste(utils::head(pair.txt, 10), collapse = "; "),
            if (length(pair.txt) > 10) {
                paste0("and ", length(pair.txt) - 10, " more")
            } else {
                ""
            },
            "\n"
        ))
    }

    # RESULTS SUMMARY [approved F11]
    if (verbose >= 3) {
        th <- fstmat[cbind(index1, index2)]
        cat(report("  Populations compared:", npops, "\n"))
        cat(report("  Pairwise comparisons:", npairs, "\n"))
        cat(report("  Weir & Cockerham theta: min",
                   round(min(th, na.rm = TRUE), 4), ", max",
                   round(max(th, na.rm = TRUE), 4), ", mean",
                   round(mean(th, na.rm = TRUE), 4), "\n"))
        cat(report("  Non-finite pairs:", length(bad.pairs), "\n"))
        if (nboots > 1) {
            cat(report("  Bootstrap replicates over loci:", nboots, "\n"))
            cat(report("  P-values are one-tailed bootstrap fractions;",
                       "no correction for multiple testing is applied\n"))
        }
    }

    # FLAG SCRIPT END

    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }

    # RETURN
    return(out)
}
