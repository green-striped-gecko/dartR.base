#' @name gl.load
#' @title Loads an object from compressed binary format produced by gl.save()
#' @family io

#' @description
#' This is a wrapper for readRDS()

#' The function loads the object from the current workspace and returns the
#'  gl object.

#' @param file Name of the file to receive the binary version of the object
#' [required].
#' @param compliance Whether to make compliance check [default FALSE]. 
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' 
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' gl.save(testset.gl,file.path(tempdir(),'testset.rds'))
#' gl <- gl.load(file.path(tempdir(),'testset.rds'))
#' @seealso \code{\link{gl.save}}
#' 
#' @export
#' @return The loaded object

gl.load <- function(file,
                    compliance = FALSE,
                    verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)
    
    x <- readRDS(file)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    cat(report("  Loaded object of type", datatype, "from", file, "\n"))

    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    invisible(x) 
    
}
