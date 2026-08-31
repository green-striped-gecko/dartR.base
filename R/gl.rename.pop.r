#' @name gl.rename.pop
#' @title Renames a population in a genlight object
#' @family data manipulation

#' @description
#' This script renames a nominated population in a genlight object.

#' @details
#' Individuals are assigned to populations based on the specimen metadata data
#' file (csv) used with gl.read.dart().
#'
#' The nominated population must exist in the genlight object, and the new
#' name must not already be in use — renaming to an existing population name
#' is a fatal error, to guard against inadvertently merging populations. To
#' amalgamate populations, use \code{\link{gl.recode.pop}}.
#'
#' The script returns a genlight object with the new population name.

#' @param x Name of the genlight object containing SNP genotypes or Tag P/A
#' data (SilicoDArT) [required].
#' @param old Name of population to be changed [required].
#' @param new New name for the population [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @return A genlight object with the new population name.

#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#'    gl <- gl.rename.pop(testset.gl, old='EmsubRopeMata', new='Outgroup')

#' @seealso \code{\link{gl.recode.pop}} for bulk recoding or amalgamating
#' populations

#' @export

gl.rename.pop <- function(x,
                         old = NULL,
                         new = NULL,
                         verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)

    # SCRIPT SPECIFIC ERROR TESTING

    if (is.null(old)) {
        stop(error("Fatal Error: A population to be renamed must be specified\n"))
    }
    if (is.null(new)) {
        stop(error("Fatal Error: A new population label must be specified\n"))
    }
    if (is.null(pop(x)) |
        is.na(length(pop(x))) | length(pop(x)) <= 0) {
        stop(error("Fatal Error: Population names not detected\n"))
    }
    if (!(old %in% popNames(x))) {
        stop(error(
            "Fatal Error: Population", old, "is not present in the dataset.
            Populations present are:",
            paste(popNames(x), collapse = ", "), "\n"
        ))
    }
    if (new %in% popNames(x)) {
        stop(error(
            "Fatal Error: Population", new, "already exists in the dataset.
            Renaming to an existing population name would merge the two
            populations; use gl.recode.pop() to amalgamate populations\n"
        ))
    }

    if (verbose >= 2) {
            cat(report("  Renaming", old, "as", new, "\n"))
    }

    # DO THE JOB

    levels(pop(x))[levels(pop(x)) == old] <- new

    # ADD TO HISTORY
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()

    # FLAG SCRIPT END

    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }

    return(invisible(x))
}
