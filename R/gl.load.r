#' @name gl.load
#' @title Loads an object from compressed binary format produced by gl.save()
#' @family io

#' @description
#' This is a wrapper for readRDS()

#' The function loads the object from the current workspace, checks if it is a
#' dartR genlight object, converts it if it is not, and returns the
#' gl object. A compliance check can be requested.

#' @param file Name of the file to receive data
#' [required].
#' @param compliance Whether to undertake a compliance check [default FALSE]. 
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' 
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
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
                     build = "v.2024.1",
                     verbose = verbose)
    
    x <- readRDS(file)
    
    # CHECK DATATYPE
    
    if (!is(x, "dartR")) {
      class(x) <- "dartR"  
      if (verbose>2) {
        cat(warn("Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"))
        cat(warn("                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR, 
                 please use function dartR2gl\n"))
      }
    }
    
    datatype <- utils.check.datatype(x, verbose = verbose)
    cat(report("  Loaded object of type", datatype, "from", file, "\n"))

    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    invisible(x) 
    
}
