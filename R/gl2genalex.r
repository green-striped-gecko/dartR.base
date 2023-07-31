#' @name gl2genalex
#' @title Converts a genlight object into a format suitable for input to genalex
#' @family linker

#' @description
#' The output csv file contains the snp data and other relevant lines suitable
#'  for genalex. This function is a wrapper for  \link[poppr]{genind2genalex}
#'  (package poppr).

#' @references
#' Peakall, R. and Smouse P.E. (2012) GenAlEx 6.5: genetic analysis
#' in Excel. Population genetic software for teaching and research-an update.
#' Bioinformatics 28, 2537-2539.
#' http://bioinformatics.oxfordjournals.org/content/28/19/2537

#' @param x Name of the genlight object containing SNP data [required].
#' @param outfile Name of the output file (including extension)
#' [default 'genalex.csv'].
#' @param outpath Path where to save the output file [default global working 
#' directory or if not specified, tempdir()].
#' @param overwrite If FALSE and filename exists, then the file will not be
#' overwritten [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end;
#' 2, progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' 
#' @author Custodian: Luis Mijangos, Author: Katrin Hohwieler, wrapper Arthur
#' Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#' \donttest{
#' gl2genalex(testset.gl, outfile='testset.csv', outpath=tempdir())
#' }
#' 
#' @export
#' @return  returns no value (i.e. NULL)

gl2genalex <- function(x,
                       outfile = "genalex.csv",
                       outpath = NULL,
                       overwrite = FALSE,
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
    
    # CHECK IF PACKAGES ARE INSTALLED
    pkg <- "poppr"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    # DO THE JOB
    
    gind <- gl2gi(x, verbose = 0)
    poppr::genind2genalex(
        gind,
        filename = outfilespec,
        sequence = TRUE,
        overwrite = overwrite
    )
    
    if (verbose > 2) {
        cat(report(paste(
            "    Records written to", outfile, ":", nInd(x), "\n"
        )))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(NULL)
    
}
