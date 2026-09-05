#' @name gl2geno
#' @title Converts a genlight object to geno format from package LEA
#' @family linker

#' @description
#' The function converts a genlight object (SNP or presence/absence
#'  i.e. SilicoDArT data) into a file in the 'geno' and the 'lfmm' formats from 
#'  (package LEA).
#'  
#' @details
#' Loci scored as missing (NA) across all individuals are removed before
#' writing, with a warning at verbose >= 1. Missing data are coded 9 in
#' both output files.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#'  (SilicoDArT) data [required].
#' @param outfile File name of the output file [default 'gl_geno'].
#' @param outpath Path where to save the output file [default global working
#' directory or if not specified, tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @return  returns no value (i.e. NULL)
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # SNP data
#' gl2geno(testset.gl, outpath=tempdir())
#' # Tag P/A data
#' gl2geno(testset.gs, outpath=tempdir())
#'
#' @export

gl2geno <- function(x,
                    outfile = "gl_geno",
                    outpath = NULL,
                    verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # SET WORKING DIRECTORY
    outpath <- gl.check.wd(outpath,verbose=0)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    hold <- x
    
    x <- gl.filter.allna(x, verbose = 0)
    
    if (nLoc(hold) > nLoc(x) & verbose > 0) {
        cat(warn(
            "  Loci with all missing data has been removed for conversion. \n"
        ))
    }
    
    # DO THE JOB
    
    outfilespec <- file.path(outpath, outfile)

    genotype <- as.matrix(x)
    genotype[is.na(genotype)] <- 9    
    # write lfmm
    write.table(
        genotype,
        paste(outfilespec, ".lfmm", sep = ""),
        col.names = FALSE,
        row.names = FALSE,
        sep = " "
    )
    # write geno
    write.table(
        t(genotype),
        paste(outfilespec, ".geno", sep = ""),
        col.names = FALSE,
        row.names = FALSE,
        sep = ""
    )
    
    if (verbose > 0) {
        cat(report(
            "  Output files:",
            paste0(outfilespec, ".geno"),
            paste0(outfilespec, ".lfmm"),
            "\n"
        ))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    invisible(NULL)
}
