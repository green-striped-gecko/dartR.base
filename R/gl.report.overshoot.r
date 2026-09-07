#' @name gl.report.overshoot
#' @title Reports loci for which the SNP has been trimmed from the sequence tag
#'  along with the adaptor
#' @family matched report

#' @description
#' This function checks the position of the SNP within the trimmed sequence tag
#' and identifies those for which the SNP position is outside the trimmed
#' sequence tag. This can happen, rarely, when the sequence containing the SNP
#' resembles the adaptor.

#' @details
#' The SNP genotype can still be used in most analyses, but functions like
#' gl2fasta() will present challenges if the SNP has been trimmed from the
#' sequence tag.

#' @param x Name of the genlight object [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @return An unaltered genlight object

#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' gl.report.overshoot(testset.gl)

#' @seealso \code{\link{gl.filter.overshoot}}

#' @export

gl.report.overshoot <- function(x,verbose = NULL) {

    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, accept = c("genlight", "SNP"), verbose = verbose)

    if (length(x@other$loc.metrics$TrimmedSequence) != nLoc(x)) {
        stop(
            error(
                "Fatal Error: Data must include Trimmed Sequences for each loci
                in a column called 'TrimmedSequence'
               in the @other$loc.metrics slot.\n"
            )
        )
    }

    if (length(x@other$loc.metrics$SnpPosition) != nLoc(x)) {
        stop(error(
            "Fatal Error: Data must include position information for each
            loci.\n"
        ))
    }

    # DO THE JOB

    if (verbose >= 2) {
        cat(report(
            "  Identifying loci for which the SNP has been trimmed with the
            adaptor\n"
        ))
    }

    trimmed <- as.character(x@other$loc.metrics$TrimmedSequence)
    snpos <- x@other$loc.metrics$SnpPosition
    # Shift the index for snppos to start from 1 not zero
    snpos <- snpos + 1
    # Pull those loci for which the SNP position is greater than the tag length
    os <- which(snpos > nchar(trimmed))

    # PRINTING OUTPUTS Report the number of such loci
    if (verbose >= 1) {
        if (length(os) > 0) {
            cat("  No. of loci with SNP falling outside the trimmed sequence:",
                length(os),
                "\n")
            cat(paste(locNames(x)[os], collapse = ", "))
            cat("\n")
        } else {
            cat("  There were no loci with SNP falling outside the trimmed sequence\n")
        }
    }

    # FLAG SCRIPT END
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n\n"))
    }

    # RETURN
    invisible(x)
}
