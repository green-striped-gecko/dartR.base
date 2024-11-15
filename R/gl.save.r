#' @name gl.save
#' @title Saves an object in compressed binary format for later rapid retrieval
#' @family io

#' @description
#' This is a wrapper for saveRDS().

#' The script saves the object in binary form to the current workspace and
#' returns the input gl object.

#' @param x Name of the genlight object containing SNP genotypes [required].
#' @param file Name of the file to receive the binary version of the object
#' [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' 
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' gl.save(testset.gl,file.path(tempdir(),'testset.rds'))
#' @seealso \code{\link{gl.load}}
#' 
#' @export
#' @return The input object

gl.save <- function(x,
                    file,
                    verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2024.1",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # DO THE JOB
    
    attributes(class(x))<- list(package="dartR.base")
    
    saveRDS(x, file)
    cat(report("  Saved object of type", datatype, "to a compressed RDA file\n"))
    cat(report("  Load again using function gl.load\n"))
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    invisible(x)
    
}
