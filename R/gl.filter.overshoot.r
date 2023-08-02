#' @name gl.filter.overshoot
#' @title Filters loci for which the SNP has been trimmed from the sequence tag
#' along with the adaptor
#' @family matched filter

#' @description
#' This function checks the position of the SNP within the trimmed sequence tag
#' and identifies those for which the SNP position is outside the trimmed
#'  sequence tag. This can happen, rarely, when the sequence containing the SNP
#'  resembles the adaptor.

#' The SNP genotype can still be used in most analyses, but functions like
#' gl2fasta() will present challenges if the SNP has been trimmed from the
#' sequence tag.

#' Not fatal, but should apply this filter before gl.filter.secondaries, for
#' obvious reasons.

#' @param x Name of the genlight object [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].

#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' result <- gl.filter.overshoot(testset.gl, verbose=3)
#' 
#' @export
#' @return A new genlight object with the recalcitrant loci deleted

gl.filter.overshoot <- function(x,
                                verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, accept = c("genlight", "SNP"), verbose = verbose)
    
    # STANDARD ERROR CHECKING
    
    # if (datatype == "SilicoDArT") {
    #     stop(
    #         error(
    #             "  Detected Presence/Absence (SilicoDArT) data. Please supply a 
    #             SNP dataset\n"
    #         )
    #     )
    # }
    
    # SCRIPT SPECIFIC ERROR CHECKING
    
    if (length(x@other$loc.metrics$TrimmedSequence) != nLoc(x)) {
        stop(
            error(
                "Fatal Error: Data must include Trimmed Sequences for each loci 
                in a column called 'TrimmedSequence' in the @other$loc.metrics 
                slot.\n"
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
    x2 <- x[, snpos > nchar(trimmed)]
    # Report the number of such loci
    if (verbose >= 3) {
        cat("  No. of loci with SNP falling outside the trimmed sequence:",
            nLoc(x2),
            "\n")
        if (nLoc(x2) > 0) {
            cat(paste0(locNames(x2), ","))
            cat("\n")
        }
    }
    
    if (verbose >= 2) {
        cat(report("  Deleting:\n"))
        cat(locNames(x2))
    }
    # extracting indexes of loci to keep
    index <- which((snpos <= nchar(trimmed)) == TRUE)
    
      # loci to keep
      x2 <- x[, index]
      # updating loc.metrics
      x2@other$loc.metrics <- x@other$loc.metrics[index, ]
    
    # # SAVE INTERMEDIATES TO TEMPDIR
    # if (save2tmp) {
    #     match_call <-
    #         paste0(names(match.call()),
    #                "_",
    #                as.character(match.call()),
    #                collapse = "_")
    #     
    #     temp_table <-
    #         tempfile(pattern = paste0(
    #             "Table",
    #             paste0(
    #                 names(match.call()),
    #                 "_",
    #                 as.character(match.call()),
    #                 collapse = "_"
    #             ),
    #             "_"
    #         ))
    #     saveRDS(data.frame(locNames = locNames(x2)), file = temp_table)
    #     if (verbose >= 2) {
    #         cat(report(
    #             "  Saving the overshot loci to the current session tempfile\n"
    #         ))
    #         cat(
    #             report(
    #                 "  NOTE: Retrieve output files from tempdir using 
    #                 gl.list.reports() and gl.print.reports()\n"
    #             )
    #         )
    #     }
    # }
    
    # ADD TO HISTORY
    nh <- length(x2@other$history)
    x2@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x2)
}
