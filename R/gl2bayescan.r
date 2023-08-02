#' @name gl2bayescan
#' @title Converts a genlight object into a format suitable for input to Bayescan
#' @family linker
#' 
#' @description
#' The output text file contains the SNP data and relevant BAyescan command
#'  lines to guide input.
#'  
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file (including extension)
#' [default bayescan.txt].
#' @param outpath Path where to save the output file [default global working 
#' directory or if not specified, tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' 
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#' out <- gl2bayescan(testset.gl, outpath = tempdir())
#' 
#' @references
#' Foll M and OE Gaggiotti (2008) A genome scan method to identify selected loci
#'  appropriate for both dominant and codominant markers: A Bayesian
#'   perspective. Genetics 180: 977-993.
#'   
#' @export
#' @return returns no value (i.e. NULL)

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
    datatype <- utils.check.datatype(x, verbose = verbose)
    
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
    
    # Create the bayescan input file
    if (verbose >= 2) {
        cat(report(
            paste("Writing text input file for Bayescan", outfilespec, "\n")
        ))
    }
    sink(outfilespec)
    
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
    
    return(NULL)
    
}
