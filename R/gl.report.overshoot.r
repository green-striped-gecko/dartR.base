#' @name gl.report.overshoot
#' @title Reports loci for which the SNP has been trimmed from the sequence tag
#'  along with the adaptor
#' @family matched report
#' 
#' @description
#' This function checks the position of the SNP within the trimmed sequence tag
#' and identifies those for which the SNP position is outside the trimmed
#' sequence tag. This can happen, rarely, when the sequence containing the SNP
#' resembles the adaptor.
#' 
#' @param x Name of the genlight object [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#' @details
#' The SNP genotype can still be used in most analyses, but functions like
#' gl2fasta() will present challenges if the SNP has been trimmed from the 
#' sequence tag.

#' Resultant ggplot(s) and the tabulation(s) are saved to the session's
#' temporary directory.

#' @author Custodian: Arthur Georges -- Post to 
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' gl.report.overshoot(testset.gl)
#' 
#' @seealso \code{\link{gl.filter.overshoot}}

#' @export
#' @return An unaltered genlight object

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
    xx <- x[, snpos > nchar(trimmed)]
    
    # PRINTING OUTPUTS Report the number of such loci
    cat("  No. of loci with SNP falling outside the trimmed sequence:",
        nLoc(xx),
        "\n")
    if (nLoc(xx) > 0) {
        cat(paste0(locNames(xx), sep = ","))
        cat("\n")
    }
    
    df <- data.frame(locNames = locNames(xx))
    
    # FLAG SCRIPT END
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n\n"))
    }
    
    # RETURN
    invisible(x)
}
