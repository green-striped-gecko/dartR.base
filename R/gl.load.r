#' @name gl.load
#' @title Loads an object from compressed binary format produced by gl.save()
#' @family io

#' @description
#' This is a wrapper for readRDS()

#' The function reads the object from the nominated file, checks if it is a
#' dartR genlight object, converts it if it is not, and returns the
#' gl object. A compliance check can be requested.
#'
#' Objects of class dartR saved by the legacy dartR package are converted to
#' dartR.base objects without loading the legacy package (which would mask
#' the dartR.base functions), and also when it is not installed. A
#' file-backed (FBM) object saved with saveRDS() instead of gl.save() can
#' only be loaded while its backing file (.bk) exists; otherwise gl.load()
#' stops with an error.

#' @param file Name of the file from which to read the saved object
#' [required].
#' @param fbm Whether to convert dartR-gen object to dartR-fbm object [default
#' FALSE].
#' @param compliance Whether to undertake a compliance check [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
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
    funname <- utils.flag.start(func = funname,
                     build = "v.2024.1",
                     verbose = verbose)
    
    if (!file.exists(file)) {
        stop(error("Fatal Error: file", file, "not found\n"))
    }

    x <- readRDS(file)

    # An S4 object records the package that defined its class. Objects of
    # class "dartR" saved by the legacy dartR package (or any package other
    # than dartR.base) point to that package, and the first is() on them
    # makes R load and attach it, which masks the dartR.base functions for
    # the rest of the session (or fails if it is not installed). The class
    # is re-pointed to dartR.base before any class check; .as_dartR() below
    # then adds the slots the legacy class lacks (@fbm).
    legacy <- FALSE
    cl <- attr(x, "class", exact = TRUE)
    if (isS4(x) && identical(as.character(cl), "dartR") &&
        !identical(attr(cl, "package"), "dartR.base")) {
      legacy <- TRUE
      attr(x, "class") <- structure("dartR", package = "dartR.base")
    }

    # An FBM-backed object saved with saveRDS() instead of gl.save() keeps
    # only a reference to its backing file (.bk), usually a temporary file
    # that is gone in a new session; gl.save() stores the genotypes instead
    fbm_in <- if (isS4(x) && methods::.hasSlot(x, "fbm")) x@fbm else NULL
    if (!is.null(fbm_in)) {
      bk <- tryCatch(fbm_in$backingfile, error = function(e) NA_character_)
      if (is.na(bk) || !file.exists(bk)) {
        stop(error("Fatal Error: file", file, "holds a file-backed (FBM)",
                   "object whose backing file", bk, "no longer exists, so",
                   "its genotypes are lost. It was most likely saved with",
                   "saveRDS(); save FBM objects with gl.save(), which stores",
                   "the genotypes\n"))
      }
    }

    if (!is(x, "genlight")) {
        stop(error(
            "Fatal Error: file", file,
            "does not contain a genlight object\n"
        ))
    }

    # CHECK DATATYPE

    if (legacy && verbose >= 2) {
      cat(warn("  Warning: Object of class dartR saved by another package",
               "(e.g. the legacy dartR) encountered. Converted to a",
               "dartR.base object\n"))
    }
    if (!is(x, "dartR")) {
      if (verbose>=2) {
        cat(warn("Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"))
        cat(warn("                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR,
                 please use function dartR2gl\n"))
      }
    }
    # Also repairs dartR objects saved without the fbm slot (made with
    # class<- or by older versions)
    x <- .as_dartR(x)

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
      x <- gl.gen2fbm(x, verbose = 0)
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
