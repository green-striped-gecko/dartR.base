#' @name gl2phylip
#' @title Creates a Phylip input distance matrix from a genlight (SNP)
#'  object
#' @family linkers
#'
#' @description
#' This function calculates and returns a matrix of Euclidean distances between
#' populations and produces an input file for the phylogenetic program Phylip
#' (Joe Felsenstein).

#' @details
#' Distances are computed from allele frequencies, locus by population: the
#' mean allele dosage of the scored individuals divided by 2. Missing
#' genotypes are excluded from each frequency; if a population has no scored
#' individuals at a locus, the cell is left NA and the affected distances
#' rest on the remaining loci (a warning is given).
#'
#' Population names are truncated to the 10-character phylip name field; a
#' warning is given when truncation renders two populations identical.
#'
#' If bstrap > 1, the loci are resampled with replacement for each bootstrap
#' replicate and the replicate distance matrices are appended to the output
#' file. The returned matrix is always the distance matrix computed from the
#' observed data.

#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile Name of the file to become the input file for phylip
#'  [default "phyinput.txt"].
#' @param outpath Path where to save the output file [default global working
#' directory or if not specified, tempdir()].
#' @param bstrap Number of bootstrap replicates [default 1].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#' @return Matrix of Euclidean distances between populations, computed from
#' the observed data.
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \donttest{
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' result <- gl2phylip(testset.gl, outfile='test.txt', outpath=tempdir(), bstrap=10)
#' }
#'
#' @import utils
#' @importFrom stats dist
#' @export

gl2phylip <- function(x,
                      outfile = "phyinput.txt",
                      outpath = NULL,
                      bstrap = 1,
                      verbose = NULL) {

    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # SET WORKING DIRECTORY (F5)
    outpath <- gl.check.wd(outpath, verbose = 0)
    outfilespec <- file.path(outpath, outfile)

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)

    # CHECK DATATYPE
    # SNP only: the frequency computation assumes diploid dosage (F4)
    datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

    # STANDARD ERROR CHECKING

    # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) |
        is.na(length(pop(x))) | length(pop(x)) <= 0) {
        if (verbose >= 2) {
            cat(
                warn(
                    "  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n"
                )
            )
        }
        pop(x) <- array("pop1", dim = nInd(x))
        pop(x) <- as.factor(pop(x))
    }

    # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose = 0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {
        cat(warn("  Warning: genlight object contains monomorphic loci\n"))
    }

    # DO THE JOB

    # Convert gl object to a matrix of allele frequencies, locus by population
    if (verbose >= 2) {
        cat(report(
            "Converting to a matrix of frequencies, locus by populations\n"
        ))
    }
    # Densify once; the bootstrap loop resamples these columns (F2, STY2)
    mat <- as.matrix(x)
    # Missing genotypes are excluded from each frequency (F3)
    t <- apply(mat, 2, tapply, pop(x), function(e)
        mean(e, na.rm = TRUE) / 2)
    # A population with no scored individuals at a locus yields NaN; the
    # affected distances rest on the remaining loci (F3, VRB4)
    empty.cells <- sum(is.nan(t))
    if (empty.cells > 0 && verbose >= 1) {
        cat(warn(
            "  Warning:", empty.cells,
            "population x locus combinations have no scored genotypes;",
            "the distances concerned are computed from the remaining loci\n"
        ))
    }
    # Compute Euclidean distance
    if (verbose >= 2) {
        cat(report("Computing Euclidean distances\n"))
    }
    d <- round(as.matrix(dist(t)), 4)
    row.names(d) <- c(paste(row.names(d), "          "))
    row.names(d) <- substr(row.names(d), 1, 10)
    # Truncation to the 10-character phylip field can collide names (F7)
    if (anyDuplicated(row.names(d)) && verbose >= 1) {
        cat(warn(
            "  Warning: population names are not unique after truncation",
            "to the 10-character phylip name field:",
            paste(unique(row.names(d)[duplicated(row.names(d))]),
                  collapse = ", "), "\n"
        ))
    }

    # Output phylip data file
    if (verbose >= 2) {
        cat(report("Writing the Phylip input file", outfilespec, "\n"))
        if (bstrap > 1) {
            cat(report(
                "Repeating calculations for",
                bstrap,
                "iterations\n"
            ))
        }
    }
    npops <- length(levels(factor(pop(x))))
    # Release the sink if anything below fails, without disturbing any
    # sink already open in the caller, e.g. capture.output (F6)
    sink.base <- sink.number()
    sink(outfilespec)
    on.exit(while (sink.number() > sink.base) sink(NULL), add = TRUE)
    cat(c("   ", npops, "\n"))
    for (i in 1:npops) {
        cat(row.names(d)[i], d[i,], "\n")
    }

    # Check if bootstrap replicates are required
    if (bstrap > 1) {
        # Repeat for each bootstrap replicate
        for (j in (2:bstrap)) {
            # Resample the loci with replacement at the matrix level:
            # genlight column subsetting with duplicated indices corrupts
            # missing genotypes (adegenet SNPbin defect) (F2)
            idx <- sample.int(nLoc(x), size = nLoc(x), replace = TRUE)

            # Convert to a matrix of allele frequencies, locus by population
            t.rep <- apply(mat[, idx, drop = FALSE], 2, tapply, pop(x),
                           function(e)
                               mean(e, na.rm = TRUE) / 2)

            # Compute Euclidean distance
            d.rep <- round(as.matrix(dist(t.rep)), 4)
            row.names(d.rep) <- c(paste(row.names(d.rep), "          "))
            row.names(d.rep) <- substr(row.names(d.rep), 1, 10)

            # Output phylip data file
            cat(c("   ", npops, "\n"))
            for (i in 1:npops) {
                cat(row.names(d.rep)[i], d.rep[i,], "\n")
            }
        }
    }
    sink()
    if (verbose >= 2) {
        cat(report("Closing output file", outfile, "\n"))
    }

    # FLAG SCRIPT END

    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }

    # Return the observed-data distance matrix, not a bootstrap replicate (F1)
    return(d)
}
