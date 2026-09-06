#' @name gl.amova
#' @title Performs AMOVA using genlight data
#' @family basic statistics

#' @description
#' This script performs an AMOVA based on the genetic distance matrix from
#' stamppNeisD() [package StAMPP] using the amova() function from the package
#' PEGAS for exploring within and between population variation. For detailed
#' information use their help pages: ?pegas::amova, ?StAMPP::stamppAmova. Be
#' aware due to a conflict of the amova functions from various packages I had
#' to 'hack' StAMPP::stamppAmova to avoid a namespace conflict.

#' @details
#' The analysis is one-level only: the populations assigned in the pop slot
#' are the single stratum. There is no provision for a nested hierarchy
#' (for example regions containing populations). At least two populations
#' are required.
#'
#' A distance matrix supplied via the distance argument must cover exactly
#' the individuals in x and must be labelled with the individual names
#' (dimnames for a matrix, Labels for a dist object). The matrix is aligned
#' to indNames(x) by those labels before use, so a matrix whose rows are
#' ordered differently from x is reordered correctly; a matrix without
#' labels, or with labels that do not match indNames(x), is a fatal error.
#'
#' Significance is assessed by the one-tailed permutation test of
#' pegas::amova(): the reported P.value is k/(permutations + 1), where k is
#' the number of permuted variance components at least as large as the
#' observed one, so a P.value of 0 is attainable. There is no seed argument;
#' for a reproducible test, call set.seed() before calling this function.
#' Set permutations = 0 to skip the test and return the variance components
#' alone.

#' @param x Name of the genlight containing the SNP genotypes, with
#' population information (at least two populations) [required].
#' @param distance Distance matrix (class matrix or dist) between
#' individuals, labelled with the individual names (if not provided NeisD
#' from StAMPP::stamppNeisD is calculated) [default NULL].
#' @param permutations Number of permutations to perform for hypothesis
#' testing [default 100]. Please note should be set to 1000 for analysis.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @return An object of class 'amova' (package pegas): a list with a table
#' (tab) of sums of square deviations (SSD), mean square deviations (MSD)
#' and degrees of freedom; the variance components (varcomp; a data.frame
#' of sigma2 and P.value when permutations > 0, otherwise a vector of
#' sigma2); the variance coefficients (varcoef); and the matched call.
#' Phi statistics are displayed when the object is printed.
#'
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' #permutations should be higher, here set to 1 because of speed
#' if (isTRUE(getOption("dartR_fbm"))) bandicoot.gl <- gl.gen2fbm(bandicoot.gl)
#' out <- gl.amova(bandicoot.gl, permutations=1)
#'
#' @export


gl.amova <- function(x,
                     distance = NULL,
                     permutations = 100,
                     verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x,
                                     accept = c("SNP", "SilicoDArT"),
                                     verbose = verbose)

    # Reject SilicoDArT before any work is done; Nei's D from stamppNeisD
    # assumes 0/1/2 SNP dosages and has no validated basis for
    # presence-absence data
    if (datatype == "SilicoDArT") {
        stop(
            error(
                "Fatal Error: Function gl.amova works only with SNP data; Nei's genetic distance is not defined for SilicoDArT presence-absence data\n"
            )
        )
    }

    # CHECK IF PACKAGES ARE INSTALLED
    pkg <- "pegas"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
        stop(error(
            "Package",
            pkg,
            " needed for this function to work. Please install it.\n"
        ))
    }

    # FUNCTION SPECIFIC ERROR CHECKING

    # AMOVA partitions variance among and within populations, so at least
    # two populations are required; a single population yields an all-NaN
    # table
    if (is.null(pop(x)) || nPop(x) < 2) {
        stop(
            error(
                "Fatal Error: AMOVA requires at least two populations; this object has",
                nPop(x),
                "\n  Assign populations to the genlight object before calling gl.amova\n"
            )
        )
    }

    #!# intermediate fbm fix
    if (!is.null(.fbm_or_null(x))) x <- gl.fbm2gen(x)

    if (is.null(distance)) {
        class(x) <- "genlight"  #needs to be genlight due to stampp
        dd <- StAMPP::stamppNeisD(x, FALSE)
        if (is.null(rownames(dd))) rownames(dd) <- indNames(x)
        if (is.null(colnames(dd))) colnames(dd) <- rownames(dd)
    } else {
        # Validate the supplied distance: it must cover exactly the
        # individuals in x and carry their names as labels, and it is
        # aligned to indNames(x) by those labels -- a mis-ordered matrix
        # applied positionally would silently misassign every distance
        if (inherits(distance, "dist")) {
            n.d <- attr(distance, "Size")
            dd <- as.matrix(distance)
        } else if (is.matrix(distance)) {
            if (nrow(distance) != ncol(distance)) {
                stop(
                    error(
                        "Fatal Error: the supplied distance matrix is not square:",
                        nrow(distance), "x", ncol(distance), "\n"
                    )
                )
            }
            n.d <- nrow(distance)
            dd <- distance
        } else {
            stop(
                error(
                    "Fatal Error: distance must be a dist object or a square matrix of distances between individuals; found",
                    class(distance)[1], "\n"
                )
            )
        }
        if (n.d != nInd(x)) {
            stop(
                error(
                    "Fatal Error: the supplied distance covers",
                    n.d,
                    "individuals; the genlight object has",
                    nInd(x), "\n"
                )
            )
        }
        labels.d <- rownames(dd)
        if (is.null(labels.d)) labels.d <- colnames(dd)
        if (is.null(labels.d)) {
            stop(
                error(
                    "Fatal Error: the supplied distance carries no individual names (dimnames for a matrix, Labels for a dist object), so it cannot be matched to the genlight object. Label it with indNames(x)\n"
                )
            )
        }
        if (anyDuplicated(labels.d) || anyDuplicated(indNames(x)) ||
            !setequal(labels.d, indNames(x))) {
            stop(
                error(
                    "Fatal Error: the labels on the supplied distance do not match indNames(x) one-to-one; the distance cannot be aligned to the genlight object\n"
                )
            )
        }
        # Align rows and columns to the order of individuals in x
        rownames(dd) <- labels.d
        colnames(dd) <- labels.d
        dd <- dd[indNames(x), indNames(x)]
    }

    # Non-finite distances (typically from individuals with all-missing
    # genotypes) propagate to an all-NaN AMOVA table -- stop and name the
    # individuals involved
    nonfinite <- !is.finite(dd)
    if (any(nonfinite)) {
        offending <- rownames(dd)[apply(nonfinite, 1, any)]
        stop(
            error(
                "Fatal Error: the distance matrix contains missing or non-finite values involving individual(s):",
                paste(offending, collapse = ", "),
                "\n  AMOVA cannot partition variance from missing distances. Consider filtering with gl.filter.callrate(method = 'ind') and gl.filter.allna() before rerunning\n"
            )
        )
    }

    pop.names <- factor(as.character(pop(x)))
    temp <- new.env()
    assign("distance", dd, envir = temp)
    assign("pop.names", pop.names, envir = temp)
    assign("permutations", permutations, envir = temp)
    res <-
        with(temp,
             pegas::amova(distance ~ pop.names, nperm = permutations))
    rm(pop.names, temp)

    # Results summary
    if (verbose >= 3) {
        cat(report("  AMOVA (one level: among and within populations)\n"))
        print(res$tab)
        if (is.data.frame(res$varcomp)) {
            sig2 <- res$varcomp$sigma2
        } else {
            sig2 <- as.numeric(res$varcomp)
        }
        cat(report(
            "  Variance components (sigma2): among =",
            signif(sig2[1], 6),
            "; within =",
            signif(sig2[2], 6),
            "\n"
        ))
        cat(report("  Phi_ST:", signif(sig2[1] / sum(sig2), 6), "\n"))
        if (is.data.frame(res$varcomp) && "P.value" %in% names(res$varcomp)) {
            cat(report(
                "  P.value (one-tailed,",
                permutations,
                "permutations):",
                res$varcomp$P.value[1],
                "\n"
            ))
        }
    }

    # FLAG SCRIPT END

    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }

    # RETURN
    return(res)

}
