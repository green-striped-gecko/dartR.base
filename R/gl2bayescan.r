#' @name gl2bayescan
#' @title Converts a genlight object into a format suitable for input to Bayescan
#' @family linker
#'
#' @description
#' The output text file contains the SNP data and relevant BayeScan command
#' lines to guide input.
#'
#' @details
#' The output follows BayeScan's codominant (GESTE) format: for each
#' population, one row per locus giving the number of gene copies sampled,
#' the number of alleles (2), and the counts of the alternate and reference
#' alleles. A population x locus combination with no genotyped individuals
#' is written as a sample of zero gene copies; a warning reports how many
#' such rows were written. Filtering on call rate
#' (\code{gl.filter.callrate}) before export is recommended to avoid empty
#' samples.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file (including extension)
#' [default bayescan.txt].
#' @param outpath Path where to save the output file [default global working
#' directory or if not specified, tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @return returns no value (i.e. NULL), invisibly
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @references
#' Foll M and OE Gaggiotti (2008) A genome scan method to identify selected loci
#'  appropriate for both dominant and codominant markers: A Bayesian
#'   perspective. Genetics 180: 977-993.
#'
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' out <- gl2bayescan(testset.gl, outpath = tempdir())
#'
#' @export

gl2bayescan <- function(x,
                        outfile = "bayescan.txt",
                        outpath = NULL,
                        verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # SET WORKING DIRECTORY
    outpath <- gl.check.wd(outpath,verbose=0)
    outfilespec <- file.path(outpath, outfile)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)
    
    # DO THE JOB
    
    if (verbose >= 2) {
        cat(report(
            paste(
                "Extracting SNP data and creating records for each individual\n"
            )
        ))
    }
    
    # Prepare the data
    mat <- gl.allele.freq(x, percent=TRUE, by='popxloc', verbose = verbose)
    mat <- mat[order(mat$popn),]
    
    # convert to character so it can be used in the for loop
    mat$popn <- as.character(mat$popn)
    
    # Warn on population x locus cells with no genotyped individuals -- these
    # are written as samples of zero gene copies, and BayeScan's tolerance of
    # empty samples is not established [approved F2]
    zero.cells <- sum(mat$nobs == 0)
    if (zero.cells > 0 && verbose >= 1) {
        cat(warn(
            "  Warning:",
            zero.cells,
            "population x locus combinations have no genotyped individuals",
            "and are written with zero gene copies. Consider filtering on",
            "call rate before export.\n"
        ))
    }

    # Create the bayescan input file
    if (verbose >= 2) {
        cat(report(
            paste("Writing text input file for Bayescan", outfilespec, "\n")
        ))
    }
    # Restore the console if anything fails while the sink is open; the
    # explicit sink() below removes the diversion on the normal path
    sink.depth <- sink.number()
    sink(outfilespec)
    on.exit(while (sink.number() > sink.depth) sink(), add = TRUE)
    
    cat(paste0("[loci]=", nLoc(x)), "\n\n")
    cat(paste0("[populations]=", nPop(x)), "\n\n")
    
    # creating a counter to be used in the for loop
    pop_id <- 0
    # adding one parenthesis that was missing and using pop names in the for loop for (i in 1:nPop(x)) {
    for (i in as.character(unique(pop(x)))) {
        # counter
        pop_id <- pop_id + 1
        # change the variable used in the loop cat(paste0('[pop]=', i), '\n')
        cat(paste0("[pop]=", pop_id), "\n")
        # change the variable used in the loop popi <- mat[mat$popn == mat$popn[i], ]
        popi <- mat[mat$popn == i,]
        for (j in 1:length(popi$popn)) {
            cat(j,
                (2 * popi$nobs[j]),
                2,
                popi$sum[j],
                (2 * popi$nobs[j] - popi$sum[j]),
                "\n")
        }
        cat("\n")
    }
    
    sink()
    
    if (verbose >= 3) {
        cat(report(paste(
            "Records written to", outfilespec, "\n"
        )))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(invisible(NULL))

}
