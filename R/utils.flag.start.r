#' @name utils.flag.start
#' @title A utility script to flag the start of a script
#' @family utilities

#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.

#' @param func Name of the function that is starting [required].
#' @param build Name of the build [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity]
#'  

#' @author Custodian: Arthur Georges -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#'  
#' @keywords internal
#' @export
#' @return The calling function name, invisibly
 
# Build Version v.2023.2

utils.flag.start <- function(func = NULL,
                             build = NULL,
                             verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    if (is.null(func)) {
        stop(error("Fatal Error: The calling function must be specified.\n"))
    }
    # func is usually match.call()[[1]] of the caller: a name, but a call
    # (pkg::fun) when called with the namespace, or the function itself under
    # do.call(fun, ...), which cat() cannot print
    func <- .call_name(func)
    if (verbose >= 1) {
        if (verbose == 5 && !is.null(build)) {
            cat(
                report(
                    "Starting",
                    func,
                    "\n[dartR.base vers.",
                    packageVersion("dartR.base"),
                    "Build =",
                    build,
                    "]\n"
                )
            )
        } else {
            cat(report("Starting", func, "\n"))
        }
    }
    invisible(func)
}

## Printable name of the function in the head of a call: a name, pkg::fun
## (or pkg:::fun), or a function object, as found under do.call(fun, ...),
## which is looked up by identity in its namespace or package environment.
.call_name <- function(f) {
    if (is.character(f)) return(f[1])
    if (is.name(f)) return(as.character(f))
    if (is.call(f) && length(f) == 3 &&
        as.character(f[[1]]) %in% c("::", ":::")) {
        return(as.character(f[[3]]))
    }
    if (is.function(f)) {
        env <- environment(f)
        if (!is.null(env)) {
            for (nm in ls(env, all.names = TRUE)) {
                obj <- get0(nm, envir = env, inherits = FALSE)
                if (is.function(obj) && identical(obj, f)) return(nm)
            }
        }
        return("<function>")
    }
    paste(deparse(f, width.cutoff = 500L, nlines = 1L), collapse = "")
}

## Short label for a value that is not worth storing in a call
.value_label <- function(v) {
    if (methods::is(v, "genlight")) {
        lab <- sprintf("%s: %d ind x %d loc", class(v)[1], nInd(v), nLoc(v))
    } else if (is.data.frame(v) || is.matrix(v)) {
        lab <- sprintf("%s: %d x %d", class(v)[1], NROW(v), NCOL(v))
    } else if (is.function(v)) {
        lab <- .call_name(v)
        if (!identical(lab, "<function>")) return(as.name(lab))
        lab <- "function"
    } else {
        lab <- sprintf("%s of length %d", class(v)[1], length(v))
    }
    as.name(paste0("<", lab, ">"))
}

#' @name utils.history.call
#' @title Compacts a call before it is stored in the history of an object
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' A call from match.call() holds the arguments as they were written, usually
#' names (e.g. x = gl). When the function is called through
#' do.call(fun, list(x = gl, ...)) the call holds the evaluated arguments
#' instead: the whole input object, which itself carries its history, so the
#' history grows with every such step. The returned call keeps names,
#' expressions and small values (atomic vectors of up to 20 elements), and
#' replaces any other value with a short label such as
#' `<dartR: 250 ind x 255 loc>`. A function given as the value of the call
#' head is replaced by its name.
#' @param call A call, typically match.call() of the calling function
#' [required].
#' @return The compacted call.
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.history.call <- function(call) {
    if (!is.call(call)) return(call)
    if (is.function(call[[1]])) {
        call[[1]] <- as.name(.call_name(call[[1]]))
    }
    for (i in seq_along(call)[-1]) {
        v <- call[[i]]
        keep <- is.null(v) || is.name(v) || is.call(v) ||
            is.expression(v) ||
            (is.atomic(v) && length(v) <= 20 && is.null(dim(v)))
        if (!keep) call[[i]] <- .value_label(v)
    }
    call
}
