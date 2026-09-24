#' @name gl.subsample.loc
#' @title Subsample loci from a genlight object
#' @family data manipulation
#'
#' @description
#' A function to subsample loci in a genlight object, at random (with or
#' without replacement) or by information content.
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param n Number of loci to include in the subsample [required].
#' @param replace If TRUE, sampling is with replacement; ignored when
#' \code{method = "pic"} [default TRUE]
#' @param error.check If TRUE, will undertake error checks on input parameters
#' [default TRUE]
#' @param method Method: "random", in which case the loci are sampled at
#' random; or "pic", in which case the top n loci ranked on information
#' content are chosen. Information content is AvgPIC for SNP data and PIC for
#' presence/absence (SilicoDArT) data, recalculated first if the stored
#' values are out of date [default "random"].
#' @param mono.rm Delete monomorphic loci before sampling [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#'
#' @details Retain a subset of loci at random, with or without replacement,
#' or the n most informative loci. Parameter n must be less than or equal to
#' nLoc(x) (after monomorphic loci are removed, if \code{mono.rm = TRUE}).
#'
#' Set error.check = FALSE for speedy execution in simulations
#'
#' @author Author(s): Bernd Gruber, Luis Mijangos. Custodian: Bernd Gruber
#' (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' gl2 <- gl.subsample.loc(testset.gl, n=50, replace=TRUE, verbose=3)
#' gl3 <- gl.subsample.loc(testset.gl, n=50, method="pic", verbose=3)
#' @export
#' @return Returns the subsampled genlight object

gl.subsample.loc <- function(x,
                             n,
                             replace = TRUE,
                             error.check = TRUE,
                             method = "random",
                             mono.rm = FALSE,
                             verbose = NULL) {
  if (missing(n)) {
    stop(error("  n, the number of loci to subsample, must be supplied\n"))
  }
  method <- tolower(method)
  if (!method %in% c("random", "pic")) {
    stop(error("  method must be 'random' or 'pic'\n"))
  }
  hold.history <- x@other$history

  if (error.check) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
  }

  if (mono.rm) {
    x <- gl.filter.monomorphs(x, verbose = 0)
  }

  if (error.check) {
    # FUNCTION SPECIFIC ERROR CHECKING

    if (n <= 0 | n > nLoc(x)) {
      if (verbose >= 1) {
        cat(warn("Subsample size must be in the range 1 to", nLoc(x), "\n"))
      }
      if (verbose >= 1) {
        cat(warn("  Set to", nLoc(x), "\n"))
      }
      n <- nLoc(x)
    }

    if (verbose >= 2) {
      if (method == "pic") {
        cat(report("  Subsampling the", n, "loci with the highest",
                   "information content from a", datatype, "object\n"))
      } else if (replace) {
        if (verbose >= 2) {
          cat(
            report(
              "  Subsampling",
              n,
              "loci at random from a",
              datatype,
              "object with replacement\n"
            )
          )
        }
      } else {
        if (verbose >= 2) {
          cat(
            report(
              "  Subsampling",
              n,
              "loci at random from a",
              datatype,
              "object without replacement\n"
            )
          )
        }
      }
    }
  }

  # DO THE JOB

  if (method == "pic") {
    # Rank on information content, recalculated when the stored values no
    # longer describe the individuals present
    if (!error.check) {
      datatype <- utils.check.datatype(x, verbose = 0)
    }
    pic.name <- if (datatype == "SNP") "AvgPIC" else "PIC"
    if (!isTRUE(x@other$loc.metrics.flags[[pic.name]]) ||
        is.null(x@other$loc.metrics[[pic.name]])) {
      x <- utils.recalc.avgpic(x, verbose = 0)
    }
    nums <- order(-x@other$loc.metrics[[pic.name]])[seq_len(n)]
    replace <- FALSE
  } else {
    # Subsample the genlight object
    # generate a random index value, with or without replacement
    nums <- sample(1:nLoc(x), size = n, replace = replace)
  }
  # subsample the data
  x2 <- x[, nums]
  # subsample the locus metrics [necessary because of replacement possibility]
  x2@other$loc.metrics <- x@other$loc.metrics[nums, ]

  if (replace == TRUE) {
    x2@loc.names <- make.unique(locNames(x2), sep = "_")
  }

  # Remove unused factor levels
  x2@other$loc.metrics[] <- lapply(x2@other$loc.metrics, function(x)
    if (is.factor(x))
      factor(x)
    else
      x)

  if (error.check) {
    # ADD TO HISTORY
    # Start from the input's history so that entries appended internally
    # (gl.filter.monomorphs) are not recorded
    x2@other$history <- hold.history
    nh <- length(x2@other$history)
    x2@other$history[[nh + 1]] <- match.call()

    # FLAG SCRIPT END ---------------

    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }
  }

  return(x2)
}
