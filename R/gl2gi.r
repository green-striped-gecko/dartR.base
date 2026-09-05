#' @name gl2gi
#' @title Converts a genlight object to genind object
#' @family linker

#' @description
#' Converts a genlight object containing SNP data into a genind object
#' (package adegenet), retaining individual names, populations and locus
#' metadata.

#' @details This function uses a faster version of df2genind (from the adegenet
#'  package).
#'
#'  Loci scored as missing (NA) across all individuals are removed before
#'  conversion, with a warning at verbose >= 1 -- df2genind would remove
#'  them silently, leaving the locus metadata out of register with the
#'  genotypes. When the genlight object carries no allele definitions
#'  (\code{loc.all} is NULL), placeholder alleles are fabricated (A/T for
#'  each locus, C/G for the first) and a warning is issued at
#'  verbose >= 1; the fabricated labels are placeholders, not sequence
#'  data.
#'
#' @param x A genlight object [required].
#' @param probar If TRUE, a progress bar is displayed for long loops
#' [default FALSE]. Retained for backward compatibility; the vectorised
#' conversion completes without one.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @return A genind object, with all slots filled.
#'
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) possums.gl <- gl.gen2fbm(possums.gl)
#' gl2gi(possums.gl[1:10,1:20])
#'
#' @export

gl2gi <- function(x,
                  probar = FALSE,
                  verbose = NULL) {

    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

    # FUNCTION SPECIFIC ERROR CHECKING

    # convert to genind....
    x2 <- as.matrix(x)
    if (!is.matrix(x2)) {
        # guard the single-locus/single-individual collapse (F3)
        x2 <- matrix(
            x2,
            nrow = nInd(x),
            ncol = nLoc(x),
            dimnames = list(indNames(x), locNames(x))
        )
    }

    # Remove loci scored NA across all individuals: df2genind would drop
    # them silently, leaving loc.metrics out of register with the
    # genotypes (F1)
    allna <- apply(x2, 2, function(v) all(is.na(v)))
    if (any(allna)) {
        if (verbose >= 1) {
            cat(warn(
                "  Loci with all missing data have been removed for",
                "conversion.\n"
            ))
        }
        keep <- which(!allna)
        hold <- x
        x <- x[, keep]
        # re-subset loc.metrics from the original object (DAT3)
        if (!is.null(hold@other$loc.metrics)) {
            x@other$loc.metrics <-
                hold@other$loc.metrics[keep, , drop = FALSE]
        }
        x2 <- x2[, keep, drop = FALSE]
    }

    # DO THE JOB

    if (probar) {
        pb <- txtProgressBar(
            min = 0,
            max = 1,
            style = 3,
            initial = NA
        )
    }

    if (is.null(x@loc.all)) {
        x@loc.all <- rep("A/T", nLoc(x))
        x@loc.all[1] <- "C/G"
        if (verbose >= 1) {
            cat(warn(
                "  Warning: no allele definitions in the genlight object;",
                "placeholder alleles (A/T, C/G for the first locus) have",
                "been fabricated. They are not sequence data.\n"
            ))
        }
    }

    homs1 <-
        paste(substr(x@loc.all, 1, 1), "/", substr(x@loc.all, 1, 1), sep = "")
    hets <- x@loc.all
    homs2 <-
        paste(substr(x@loc.all, 3, 3), "/", substr(x@loc.all, 3, 3), sep = "")

    # Vectorised recode (F8): index the per-locus genotype strings by
    # column in one pass rather than an elementwise double loop
    xx <- matrix("-/-", nrow = nrow(x2), ncol = ncol(x2))
    colidx <- col(x2)
    idx <- which(!is.na(x2) & x2 == 0)
    xx[idx] <- homs1[colidx[idx]]
    idx <- which(!is.na(x2) & x2 == 1)
    xx[idx] <- hets[colidx[idx]]
    idx <- which(!is.na(x2) & x2 == 2)
    xx[idx] <- homs2[colidx[idx]]

    if (probar) {
        setTxtProgressBar(pb, 1)
        close(pb)
    }

    if (verbose >= 2) {
        cat(report("Matrix converted.. Prepare genind object...\n"))
    }

    gen <-
        df2genind(
            xx[, , drop = FALSE],
            sep = "/",
            ncode = 1,
            ind.names = x@ind.names,
            pop = x@pop,
            ploidy = 2,
            NA.char = "-",
            loc.names = locNames(x)
        )  #, probar=probar)
    gen@other <- x@other

    # Keep loc.metrics in register with the loci actually present in the
    # genind (DAT3 idiom; F1)
    if (!is.null(gen@other$loc.metrics)) {
        keep <- match(locNames(gen), locNames(x))
        gen@other$loc.metrics <-
            x@other$loc.metrics[keep, , drop = FALSE]
    }

    # FLAG SCRIPT END

    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }

    # RETURN

    return(gen)

}
