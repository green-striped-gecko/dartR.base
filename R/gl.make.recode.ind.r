#' @name gl.make.recode.ind
# Preliminaries -- Parameter specifications -------------- 
#' @title Creates a proforma recode_ind file for reassigning individual
#'  (=specimen) names
#' @family data manipulation

#' @description
#' Renaming individuals may be required when there have been errors in labeling
#'  arising in the process from sample to sequencing files. There may be occasions
#'  where renaming individuals is required for preparation of figures. 
#' 
#' @param x Name of the genlight object [required].
#' @param out.recode.file File name of the output file (including extension)
#'  [default default_recode_ind.csv].
#' @param outpath Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()] 
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, 
#' progress log; 3, progress and results summary; 5, full report 
#' [default 2 or as specified using gl.set.verbosity].
#' 
#' @details
#' This function facilitates the construction of a recode table by producing a
#'  proforma file with current individual (=specimen) names in two identical
#'  columns. Edit the second column to reassign individual names. Use keyword
#'  'Delete' to delete an individual.

#'  When caution needs to be exercised because of the potential for breaking the
#'  'chain of evidence' associated with the samples, recoding individuals using
#'  a recode table (csv) can provide a clear record of the changes.

#' Use outpath=getwd() or when calling this function to direct output files 
#' to your working directory.

#' The function works with both genlight objects
#' containing SNP genotypes and Tag P/A data (SilicoDArT).

#' Apply the recoding using gl.recode.ind(). 

#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

# Examples --------------
#' @examples
#' result <- gl.make.recode.ind(testset.gl, out.recode.file ='Emmac_recode_ind.csv',outpath=tempdir())

#' @export
#' @return A vector containing the new individual names.

# Function --------------
gl.make.recode.ind <- function(x,
                               out.recode.file = "default_recode_ind.csv",
                               outpath = NULL,
                               verbose = NULL) {
  # Preliminaries -------------
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
    
    outfilespec <- file.path(outpath, out.recode.file)
    
    # DO THE JOB --------------
    
    # if (verbose >= 2) {cat(report(' Creating draft lookup table\n'))}
    mat <- cbind(indNames(x), indNames(x))
    if (verbose >= 2) {
        cat(report(
            "  Writing draft lookup table to",
            outfilespec,
            ". Edit before use\n"
        ))
    }
    write.table(
        mat,
        file = outfilespec,
        sep = ",",
        row.names = FALSE,
        col.names = FALSE
    )
    
    # FLAG SCRIPT END ---------------
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    # End block ---------------
    
    return(NULL)
}
