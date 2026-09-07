#' @name gl.load
#' @title Loads an object from compressed binary format produced by gl.save()
#' @family io

#' @description
#' This is a wrapper for readRDS()

#' The function reads the object from the nominated file, checks if it is a
#' dartR genlight object, converts it if it is not, and returns the
#' gl object. A compliance check can be requested.

#' @param file Name of the file from which to read the saved object
#' [required].
#' @param fbm Whether to convert dartR-gen object to dartR-fbm object [default
#' FALSE].
#' @param compliance Whether to undertake a compliance check [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#'
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' f <- file.path(tempdir(), 'testset.rds')
#' gl.save(testset.gl, f)
#' gl <- gl.load(f)
#'
#' @seealso \code{\link{gl.save}}
#'
#' @return The loaded object, returned invisibly
#' @export

gl.load <- function(file,
                    fbm = FALSE,
                    compliance = FALSE,
                    verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2024.1",
                     verbose = verbose)
    
    if (!file.exists(file)) {
        stop(error("Fatal Error: file", file, "not found\n"))
    }

    x <- readRDS(file)

    if (!is(x, "genlight")) {
        stop(error(
            "Fatal Error: file", file,
            "does not contain a genlight object\n"
        ))
    }

    # CHECK DATATYPE

    if (!is(x, "dartR")) {
      class(x) <- "dartR"
      if (verbose>=2) {
        cat(warn("Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"))
        cat(warn("                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR,
                 please use function dartR2gl\n"))
      }
    }

    datatype <- utils.check.datatype(x, verbose = verbose)
    if (verbose >= 2) {
        cat(report("  Loaded object of type", datatype, "from", file, "\n"))
    }

    # Undertake a compliance check only if requested
    if (compliance) {
        x <- gl.compliance.check(x, verbose = verbose)
    }

    #convert to fbm if requested
    if (fbm)    {
      x <- dartR.base::gl.gen2fbm(x)
      if (verbose >= 2) {
        cat(report("  Converted dartR-gen object to dartR-FBM object\n"))
      }
    }
    
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    invisible(x) 
    
}
