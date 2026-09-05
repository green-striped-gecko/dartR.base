#' @name gl2treemix
#' @title Converts a genlight object to a treemix input file
#' @family linker

#' @description
#' The output file contains the SNP data in the format expected by treemix --
#' see the treemix manual. The file will be gzipped before in order to be
#' recognised by treemix. Plotting functions provided with treemix will need to
#'  be sourced from the treemix download page.


#' @details
#' The output file has a header row of population names followed by one
#' row per locus, each cell holding the reference and alternate allele
#' counts for that population as 'ref,alt'.
#'
#' @param x Name of the genlight object [required].
#' @param outfile File name of the output file (including gz extension)
#' [default 'treemix_input.gz'].
#' @param outpath Path where to save the output file [default global working
#' directory or if not specified, tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @return  returns no value (i.e. NULL)
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' gl2treemix(testset.gl, outpath=tempdir())
#' 
#' @references Pickrell and Pritchard (2012). Inference of population splits and
#' mixtures from genome-wide allele frequency data. PLoS Genetics
#'  https://doi.org/10.1371/journal.pgen.1002967
#'
#' @export

gl2treemix <- function(x,
                       outfile = "treemix_input.gz",
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
    # SNP only: the allele-count arithmetic below assumes two allele
    # copies per scored individual, which is wrong for ploidy-1
    # (SilicoDArT) data
    datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)
    
    # DO THE JOB
    
    freq <- gl.allele.freq(x, percent=TRUE, by='popxloc', verbose = verbose)
    freq$ref <- freq$nobs * 2 - freq$sum
    freq$alt <- freq$sum
    freq$sum <- NULL
    freq$nobs <- NULL
    freq$nmissing <- NULL
    freq$frequency <- NULL
    freq$n <- NULL
    
    # Output the file
    
    if (verbose >= 2) {
        cat(report(
            paste(
                "    Writing results to treemix input file",
                outfilespec,
                "\n"
            )
        ))
    }
    # Restore the console and close the gz connection if the write fails
    # or is interrupted; the success path closes both explicitly below
    sink.depth <- sink.number()
    con <- gzfile(outfilespec)
    sink(con)
    on.exit({
        while (sink.number() > sink.depth) sink()
        tryCatch(close(con), error = function(e) NULL)
    }, add = TRUE)

    cat(unique(as.character(freq$popn)), "\n")
    k <- 1
    for (j in 1:nLoc(x)) {
        for (i in k:(k + nPop(x) - 1)) {
            cat(paste0(freq$ref[i], ",", freq$alt[i]), " ")
        }
        cat("\n")
        k <- k + nPop(x)
    }

    sink()
    close(con)
    if (verbose > 2) {
        cat(report(
            paste("    Records (loci) written to", outfilespec, ":",
                  nLoc(x), "for", nPop(x), "populations\n")
        ))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
        cat(important("Output file has been gzipped for input to treemix\n"))
    }
    
    invisible(NULL)

}
