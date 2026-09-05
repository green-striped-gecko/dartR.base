#' @name gi2gl
#' @title Converts a genind object into a genlight object
#' @family linker

#' @description
#' Converts a genind object (package adegenet) containing biallelic SNP
#' data into a genlight object, retaining individual names, populations
#' and the contents of the other slot.

#' @details
#' Only diploid genind objects with at most two alleles per locus can be
#' converted; loci with more than two alleles cannot be represented in a
#' genlight object and raise a fatal error naming the offending loci.
#'
#' Allele spellings are reconstructed from the genind allele names: each
#' locus receives \code{loc.all} of the form second/first allele of the
#' genind tab, so the 0/1/2 dosages and the allele pair are mutually
#' consistent. Be aware that due to ambiguity over which allele is the
#' reference, a combination of gi2gl(gl2gi(gl)) does not return an
#' identical object (loci can come back with the complementary dosage and
#' the reversed allele pair), but in terms of analysis the conversions
#' are equivalent.
#'
#' @param gi A genind object [required].
#' @param parallel Switch to deactivate parallel version. It might not be worth
#' to run it parallel most of the times [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @return A genlight object, with all slots filled.
#'
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' gi <- gl2gi(gl.filter.allna(testset.gl[1:20, 1:50], verbose = 0),
#'             verbose = 0)
#' gl <- gi2gl(gi, verbose = 0)
#'
#' @export

gi2gl <- function(gi,
                  parallel = FALSE,
                  verbose = NULL) {

    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)

    # STANDARD ERROR CHECKING

    if (!is(gi, "genind")) {
        stop(error("  Fatal Error: genind object required!\n"))
    }

    # FUNCTION SPECIFIC ERROR CHECKING

    # Only diploid genotypes can be represented as 0/1/2 dosages (F5)
    if (any(adegenet::ploidy(gi) != 2)) {
        stop(error(
            "  Fatal Error: only diploid genind objects can be converted",
            "to a genlight object.\n"
        ))
    }

    # A genlight locus is biallelic. Loci with more than two alleles
    # cannot be represented, and previously shifted every subsequent
    # locus onto the wrong tab column, silently corrupting the output
    # (F1)
    locna <- gi@loc.n.all
    if (any(locna > 2)) {
        offending <- locNames(gi)[locna > 2]
        stop(error(
            "  Fatal Error: loci with more than two alleles cannot be",
            "converted to a genlight object:",
            paste(offending, collapse = ", "),
            "\n"
        ))
    }

    # DO THE JOB

    # Index of the first tab column of each locus: step by the actual
    # allele count of the preceding locus (F1, F2; the previous walk
    # hard-coded steps of 1 or 2 and crashed on a single-locus genind)
    ccc <- cumsum(c(1, locna[-length(locna)]))

    gl <-
        new(
            "genlight",
            gi@tab[, ccc, drop = FALSE],
            pop = pop(gi),
            other = gi@other,
            ploidy = 2,
            loc.names = locNames(gi),
            ind.names = indNames(gi),
            parallel = parallel
        )

    # Reconstruct allele spellings from the genind allele names (F4).
    # The dosage counts the first-listed allele, and dartR codes 0 as
    # homozygous for the first allele of loc.all, so the pair is written
    # second/first; a locus monomorphic in the data comes back as 'a/a'.
    gl@loc.all <- unname(vapply(gi@all.names, function(a) {
        if (length(a) >= 2) {
            paste0(a[2], "/", a[1])
        } else {
            paste0(a[1], "/", a[1])
        }
    }, character(1)))

    gl <- gl.compliance.check(gl, verbose = verbose)

    if (is.null(gl@other$loc.metrics.flags$monomorphs)) {
        gl@other$loc.metrics.flags$monomorphs <- FALSE
    }

    # ADD TO HISTORY

    nh <- length(gl@other$history)
    gl@other$history[[nh + 1]] <- match.call()

    # FLAG SCRIPT END

    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }

    return(gl)
}
