#' @name gl.fdsim
#' @title Estimates the rate of false positives in a fixed difference analysis
#' @family fixed difference analysis
#'
#' @description
#' Estimates how many fixed differences between two populations would arise
#' from sampling error alone, given the observed allele frequencies and
#' sample sizes, and the probability that false positives alone reach the
#' observed count of fixed differences.
#'
#' @details
#' Loci that are missing in either population are dropped, as are loci
#' that already show a true fixed difference: minor allele frequency below
#' \code{delta} in one population and above \code{1 - delta} in the other.
#' Each remaining locus can produce a false positive.
#'
#' In each of \code{reps} replicates, a true allele frequency is drawn for
#' each population and locus by binomial sampling around the observed
#' frequency. From it, the probability that samples of the observed size
#' show a fixed difference is computed, and a false positive is drawn with
#' that probability. The count of false positives across loci is recorded.
#'
#' For allopatric populations (\code{sympatric = FALSE}) each population is
#' sampled from its own observed frequencies. For sympatric populations
#' (\code{sympatric = TRUE}) both are sampled from the pooled frequency,
#' weighted by sample size, each with its own sample size: the null
#' hypothesis that both samples come from one gene pool.
#'
#' Sample sizes are counted in alleles (twice the number of individuals)
#' for SNP data, and in individuals for SilicoDArT data.
#'
#' The p-value is the share of replicates whose count reaches the observed
#' count, \code{(sum(count >= obs) + 1) / (reps + 1)}, so the smallest
#' reportable value is \code{1 / (reps + 1)}.
#'
#' @param x Name of the genlight object containing the SNP or SilicoDArT
#' genotypes [required].
#' @param poppair Labels of two different populations for comparison, in the
#' form c(popA, popB) [required].
#' @param obs Observed number of fixed differences between the two
#' populations. If NULL, it is calculated with gl.fixed.diff [default NULL].
#' @param sympatric If TRUE, the two populations are sympatric; if FALSE,
#' allopatric [default FALSE].
#' @param reps Number of replicates in the simulation, a whole number of at
#' least 2 [default 1000].
#' @param delta The threshold value for the minor allele frequency to regard
#' the difference between two populations as fixed [default 0.02].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @return A named list of four numbers:
#' \itemize{
#'   \item observed -- the observed count of fixed differences;
#'   \item mnexpected -- the mean simulated count of false positives;
#'   \item sdexpected -- the standard deviation of the simulated count of
#'   false positives;
#'   \item prob -- the probability that false positives alone reach the
#'   observed count.
#' }
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges (Post to
#'  \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' fd <- gl.fdsim(testset.gl[, 1:100],
#'   poppair = c("EmsubRopeMata", "EmmacBurnBara"),
#'   reps = 100, verbose = 3
#' )
#'
#' @importFrom stats rbinom runif sd
#' @export

gl.fdsim <-  function(x,
                      poppair,
                      obs = NULL,
                      sympatric = FALSE,
                      reps = 1000,
                      delta = 0.02,
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

    # SCRIPT SPECIFIC CHECKS

    if (length(poppair) != 2 || poppair[1] == poppair[2]) {
        stop(error(
            "  Fatal Error: poppair must be two different population labels,",
            "e.g. poppair = c(popA, popB)\n"
        ))
    }
    for (p in poppair) {
        if (!(p %in% levels(pop(x)))) {
            stop(error("  Fatal Error: population", p,
                       "not found in the genlight object\n"))
        }
    }
    if (!is.numeric(reps) || length(reps) != 1 || is.na(reps) ||
        reps < 2 || reps != round(reps)) {
        stop(error("  Fatal Error: reps must be a whole number of at least 2\n"))
    }
    if (!is.null(obs) && (!is.numeric(obs) || length(obs) != 1 ||
                          is.na(obs) || obs < 0)) {
        stop(error("  Fatal Error: obs must be NULL or a single",
                   "non-negative number\n"))
    }

    # DO THE JOB

    # Extract the data for the two nominated populations
    pair <- gl.keep.pop(x,
                        pop.list = poppair,
                        recalc = FALSE,
                        mono.rm = TRUE,
                        verbose = 0)

    if (verbose >= 2) {
        cat(report("    Populations", poppair[1], "vs", poppair[2],
                   if (sympatric) "[sympatric]\n" else "[allopatric]\n"))
    }
    if (verbose >= 3) {
        n.ind <- table(pop(pair))[poppair]
        cat(report("    Sample sizes:",
                   paste0(poppair, " = ", n.ind, collapse = ", "), "\n"))
        cat(report("    No. of loci:", nLoc(pair), "\n"))
    }

    # Allele frequencies (as proportions) and sample sizes per locus; the
    # sampling unit is the allele for SNP data, the individual for SilicoDArT
    rf <- gl.allele.freq(pair, percent = TRUE, by = "popxloc", verbose = 0)
    rfA <- rf[rf$popn == poppair[1], ]
    rfB <- rf[rf$popn == poppair[2], ]
    pA <- rfA$frequency / 100
    pB <- rfB$frequency / 100
    ploidy <- if (datatype == "SNP") 2 else 1
    nA <- rfA$nobs * ploidy
    nB <- rfB$nobs * ploidy

    # Candidate loci: frequency known in both populations, and not already a
    # (near) fixed difference for the given delta
    keep <- !is.na(pA) & !is.na(pB) &
        !((pA < delta & (1 - pB) < delta) | (pB < delta & (1 - pA) < delta))
    pA <- pA[keep]
    pB <- pB[keep]
    nA <- nA[keep]
    nB <- nB[keep]
    nloc <- sum(keep)

    # Sympatric: both samples drawn from one gene pool
    if (sympatric) {
        pA <- pB <- (pA * nA + pB * nB) / (nA + nB)
    }

    # Calculate the observed fixed differences
    if (is.null(obs)) {
        fdmat <- gl.fixed.diff(pair, verbose = 0)
        obs <- fdmat$fd[1]
    }

    # Simulate the count of false positives
    if (verbose >= 2) {
        cat(report("  Calculating false positive rate with", reps,
                   "replications. Please be patient\n"))
    }
    falsepos <- numeric(reps)
    for (j in 1:reps) {
        # Draw a true allele frequency for each population and locus
        simA <- stats::rbinom(nloc, size = nA, prob = pA) / nA
        simB <- stats::rbinom(nloc, size = nB, prob = pB) / nB
        # Probability that samples of the observed sizes are fixed for
        # alternative alleles (Equation 5)
        pfd <- (1 - simA)^nA * simB^nB + simA^nA * (1 - simB)^nB
        # Draw whether each locus produces a false positive, and count them;
        # the spread of this count is what the p-value needs
        falsepos[j] <- sum(stats::runif(nloc) < pfd)
    }
    mn <- mean(falsepos)
    sdev <- stats::sd(falsepos)
    # Probability that false positives alone reach the observed count
    nprob <- (sum(falsepos >= obs) + 1) / (reps + 1)

    if (verbose >= 3) {
        cat(report("    Threshold minor allele frequency for a true fixed",
                   "difference:", delta, "\n"))
        cat(report("    Mean simulated count of false positives:",
                   round(mn, 2), "\n"))
        cat(report("    SD of the simulated count of false positives:",
                   round(sdev, 2), "\n"))
        cat(report("    Probability that false positives alone reach the",
                   "observed count of", obs, ":", signif(nprob, 4), "\n"))
    }

    l <- list(observed = obs,
              mnexpected = mn,
              sdexpected = sdev,
              prob = nprob)

    # FLAG SCRIPT END

    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }

    return(l)
}
