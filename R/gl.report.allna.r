#' @name gl.report.allna
#' @title Reports loci that are all NA across individuals and/or populations
#' with all NA across loci
#' @family matched report

#' @description
#' This script reports loci or individuals with all calls missing (NA),
#'  from a genlight object.
#'
#' A DArT dataset will not have loci for which the calls are scored all as
#' missing (NA) for a particular individual, but such loci can arise rarely when
#'  populations or individuals are deleted. Similarly, a DArT dataset will not
#'  have individuals for which the calls are scored all as missing (NA) across
#'  all loci, but such individuals may sneak in to the dataset when loci are
#'  deleted. Retaining individual or loci with all NAs can cause issues for
#'  several functions.
#'
#'  Also, on occasions an analysis will require that there are some loci scored
#'  in each population. Setting by.pop=TRUE will result in removal of loci when
#'  they are all missing in any one population.
#'
#' Note that loci that are missing for all individuals in a population are
#' not imputed with method 'frequency' or 'HW'. Consider
#' using the function \code{\link{gl.filter.allna}} with by.pop=TRUE.

#' @param x Name of the input genlight object [required].
#' @param by.pop If TRUE, loci that are all missing in any one population
#' are reported [default FALSE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].

#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # SNP data
#'   result <- gl.report.allna(testset.gl, verbose=3)
#' # Tag P/A data
#'   result <- gl.report.allna(testset.gs, verbose=3)

#' @family filter functions
#' @import utils patchwork
#' @export
#' @return gl.report.allna

gl.report.allna <- function(x,
                            by.pop = FALSE,
                            verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)
  
  # CHECK DATATYPE
  # datatype <- utils.check.datatype(x, accept = c("genlight", "SNP", "SilicoDArT"), verbose = verbose)
  
  # DO THE JOB
  
  if (verbose >= 2) {
    if (by.pop == FALSE) {
      cat(report(
        "  Identifying loci and individuals scored all missing (NA)\n"
      ))
    } else {
      cat(report(
        "  Identifying loci that are all missing (NA) in any one population\n"
      ))
    }
  }
  
  if (by.pop == FALSE) {
    # Consider loci
    if (verbose >= 2) {
      cat(report("  Reporting loci that are scored as all missing (NA)\n"))
    }
    na.counter <- 0
    nL <- nLoc(x)
    loc.list <- array(NA, nL)
    matrix_tmp <- as.matrix(x)
    l.names <- locNames(x)
    for (i in 1:nL) {
      row_tmp <- matrix_tmp[, i]
      if (all(is.na(row_tmp))) {
        loc.list[i] <- l.names[i]
        na.counter <- na.counter + 1
      }
    }
    if (na.counter == 0) {
        cat("  Zero loci that are missing (NA) across all individuals\n")
    } else {
      loc.list <- loc.list[!is.na(loc.list)]
        cat(
          "  ",
          na.counter,
          "loci that are missing (NA) across all individuals:",
          paste(loc.list, collapse = ", "),
          "\n"
        )
    }
    
    # Consider individuals
    if (verbose >= 2) {
      cat(report(
        "  Reporting individuals that are scored as all missing (NA)\n"
      ))
    }
    na.counter <- 0
    nI <- nInd(x)
    ind.list <- vector("list", nI)
    matrix_tmp <- as.matrix(x)
    I.names <- indNames(x)
    for (i in 1:nI) {
      row_tmp <- matrix_tmp[i, ]
      if (all(is.na(row_tmp))) {
        ind.list[i] <- I.names[i]
        na.counter <- na.counter + 1
      }
    }
    
    if (na.counter == 0) {
        cat("  Zero individuals that are missing (NA) across all loci\n")
    } else {
      ind.list <- ind.list[!is.na(ind.list)]
        cat(
          "  Individuals that are missing (NA) across all loci:",
          paste(ind.list, collapse = ", "),
          "\n"
        )
    }
  }
  
  if (by.pop == TRUE) {
    if (verbose >= 2) {
      cat(report(
        "  Reporting loci that are all missing (NA) in any one population\n"
      ))
    }
    total <- 0
    loc.list <- NULL
    for (i in 1:nPop(x)) {
      tmpop <- as.matrix(gl.keep.pop(x, popNames(x)[i], verbose = 0))
      tmpsums <- apply(tmpop, 2, function(x) {
        all(is.na(x))
      })
      tmp.list <- locNames(x)[tmpsums == TRUE]
      count <- length(tmp.list)
      cat("    ", popNames(x)[i], count, "loci with all missing data\n")
      total <- total + count
      loc.list <- c(loc.list, tmp.list)
    }
    loc.list <- unique(loc.list)
      cat("\n  Loci all NA in one or more populations:",
          length(loc.list),"\n")
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(invisible(x))
}
