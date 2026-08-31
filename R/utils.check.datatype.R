#' @name utils.check.datatype
#' @title Utility function to check the class of an object passed to a function
#' @family utilities
#' 
#' @description 
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.

#' @param x Name of the genlight object, dist matrix, data matrix, glPCA, or
#' fixed difference list (fd) [required].
#' @param accept Vector containing the classes of objects that are to be
#' accepted. The entries 'genlight' and 'dartR' accept both genotype
#' datatypes ('SNP' and 'SilicoDArT')
#' [default c('genlight','SNP','SilicoDArT','dartR')].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#' 
#' @details
#' Most functions require access to a genlight object, dist matrix, data matrix
#'  or fixed difference list (fd), and this function checks that a genlight
#'  object or one of the above has been passed, whether the genlight object is a
#'   SNP dataset or a SilicoDArT object, and reports back if verbosity is >=2.
#'   
#' This function checks the class of passed object and sets the datatype to
#' 'SNP', 'SilicoDArT', 'dist', 'matrix', 'glPca', 'fd', 'list', or
#' class(x)[1] as appropriate.

#' Note also that this function checks to see if there are individuals or loci
#' scored as all missing (NA) and if so, issues the user with a warning
#' (verbose >= 2).

#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' datatype <- utils.check.datatype(testset.gl)
#' datatype <- utils.check.datatype(as.matrix(testset.gl),accept='matrix')
#' fd <- gl.fixed.diff(testset.gl)
#' datatype <- utils.check.datatype(fd,accept='fd')
#'  
#' @return datatype (returned invisibly): 'SNP' for SNP data, 'SilicoDArT' for
#' P/A data, 'dist' for a distance matrix, 'matrix' for a data matrix, 'glPca'
#' for an ordination file, 'fd' for a fixed difference object, 'list' for a
#' list, or class(x)[1].
#' @export

utils.check.datatype <- function(x,
                                 accept = c("genlight", "SNP", "SilicoDArT","dartR"),
                                 verbose = NULL) {
    #### SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    if (is(x, "dartR") | is(x,"genlight")) {
        if (is.null(ploidy(x))) {
            stop(
                error(
                    "Fatal Error: ploidy not set in the genlight object, run gl <- gl.compliance.check(gl)\n"
                )
            )
        }
        if (verbose >= 2) {
            cat(report("  Processing genlight object"))
        }
        if (all(ploidy(x) == 1)) {
            if (verbose >= 2) {
                cat(report(" with Presence/Absence (SilicoDArT) data\n"))
            }
            datatype <- "SilicoDArT"
        } else if (all(ploidy(x) == 2)) {
            if (verbose >= 2) {
                cat(report(" with SNP data\n"))
            }
            datatype <- "SNP"
        } else {
          # Mixed or non-diploid ploidy (e.g. polyploid data) is treated
          # as SNP data, with notice
          if (verbose >= 2) {
            cat(report(" with SNP data\n"))
            cat(warn(
              "  Warning: ploidy is not uniformly 2; treating as SNP data\n"
            ))
          }
          datatype <- "SNP"
        }
        # Check for individuals or loci scoring all missing values (NA)
        # (a direct single-pass check; running a full gl.filter.allna
        # here cost ~50 ms on every function entry)
        if (verbose > 1) {
            namat <- is.na(as.matrix(x))
            if (any(colSums(!namat) == 0)) {
                cat(
                    warn(
                        "  Warning: data include loci that are scored NA across all individuals.\n  Consider filtering using gl <- gl.filter.allna(gl)\n"
                    )
                )
            }
            if (any(rowSums(!namat) == 0)) {
                cat(
                    warn(
                        "  Warning: data include individuals that are scored NA across all loci.\n  Consider filtering using gl <- gl.filter.allna(gl)\n"
                    )
                )
            }
        }
    } else if (is(x, "fd")) {
        if (is(x$gl, "genlight")) {
            # if(is.null(ploidy(x$gl))){ stop(error('Fatal Error: ploidy not set in the genlight object, run gl <-
            # gl.compliance.check(gl)\n')) }
            if (verbose >= 2) {
                cat(report("  Processing a fixed difference (fd) object"))
            }
            if (verbose >= 2) {
                if (all(ploidy(x$gl) == 1)) {
                    cat(report(
                        " with Presence/Absence (SilicoDArT) data\n"
                    ))
                } else {
                    cat(report(" with SNP data\n"))
                }
            }
        } else {
            stop(
                error(
                    "Fatal Error: Fixed Difference object expected! Check format of object\n"
                )
            )
        }
        datatype <- "fd"
    } else if (is(x, "dist")) {
        if (verbose >= 2) {
            cat(report("  Processing a distance matrix\n"))
        }
        datatype <- "dist"
    } else if (is(x, "matrix")) {
        if (verbose >= 2) {
            cat(report("  Processing a data matrix\n"))
        }
        datatype <- "matrix"
    } else if (is(x, "glPca")) {
        if (verbose >= 2) {
            cat(report("  Processing an ordination file (glPca)\n"))
        }
        datatype <- "glPca"
    } else if (is.list(x)) {
        if (verbose >= 2) {
            cat(report("  Processing a list\n"))
        }
        datatype <- "list"
    } else {
        if (verbose >= 1) {
            cat(warn("  Warning: Found object of class", class(x)[1], "\n"))
        }
        datatype <- class(x)[1]
    }

    #### CHECK WHETHER TO THROW AN ERROR ####

    # 'genlight' or 'dartR' in accept admits both genotype datatypes,
    # unless a specific genotype datatype is also listed (in which case
    # the specific listing governs, e.g. c('genlight','SNP') is SNP-only)
    if (datatype %in% c("SNP", "SilicoDArT") &&
        any(c("genlight", "dartR") %in% accept) &&
        !any(c("SNP", "SilicoDArT") %in% accept)) {
        accept <- union(accept, c("SNP", "SilicoDArT"))
    }

    if (!(datatype %in% accept)) {
        stop(
            error(
                "Fatal Error: inappropriate object passed to function, found",
                datatype,
                "expecting",
                paste(accept, collapse = " or ")
            )
        )
    }
    
    invisible(datatype)
}

