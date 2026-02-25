#' @name gl.reassign.ind
#' @title Reassign specified individuals to a (new or existing) population
#' @family data manipulation
#'
#' @description
#' Reassigns a specified set of individuals from their current populations
#' to a new population label, which may or may not already exist in the
#' genlight \{adegenet\} object.
#'
#' This function is useful when you want to manually override the population
#' membership of a subset of individuals, for example to combine individuals
#' from several small populations into a single larger one, or to pull a few
#' misassigned individuals into a corrected population.
#'
#' The function returns a genlight object with updated population assignments
#' for the specified individuals. All other individuals retain their original
#' population assignments.
#'
#' @param x Name of the genlight object containing SNP genotypes [required].
#' @param ind.list Vector specifying individuals to be reassigned. Can be a
#'  character vector of individual names (matching \code{indNames(x)}),
#'  a numeric vector of indices, or a logical vector of length \code{nInd(x)}
#'  [required].
#' @param new.pop Character string giving the population label to which the
#'  specified individuals will be reassigned. This value may or may not
#'  already exist as an entry in \code{pop(x)} [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log; 3, progress and results summary; 5, full report
#'  [default 2 or as specified using \code{gl.set.verbosity}].
#'
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#'   # SNP data
#'   popNames(testset.gl)
#'   indNames(testset.gl)[1:5]
#'   gl2 <- gl.reassign.ind(
#'     testset.gl,
#'     ind.list = indNames(testset.gl)[1:5],
#'     new.pop  = "NewGroup",
#'     verbose  = 3
#'   )
#'   popNames(gl2)
#'
#'   # Tag P/A data
#'   popNames(testset.gs)
#'   gs2 <- gl.reassign.ind(
#'     testset.gs,
#'     ind.list = 1:10,
#'     new.pop  = "Reassigned",
#'     verbose  = 3
#'   )
#'   popNames(gs2)
#' }
#'
#' @return A genlight object with the specified individuals reassigned to
#'  the new population.
#' @export
#'
gl.reassign.ind <- function(x,
                            ind.list,
                            new.pop,
                            verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2025.1",
                   verbose = verbose)
  
  # CHECK DATATYPE (for consistency with other dartR functions)
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # ---- CHECKS ON ARGUMENTS ----
  
  if (missing(ind.list) || is.null(ind.list)) {
    stop("Fatal error: ind.list must be provided\n")
  }
  
  if (missing(new.pop) || is.null(new.pop) || new.pop == "") {
    stop("Fatal error: new.pop must be a non-empty character string\n")
  }
  
  # Standardise individual indices
  n_ind <- nInd(x)
  all_inds <- indNames(x)
  
  if (is.logical(ind.list)) {
    if (length(ind.list) != n_ind) {
      stop("Fatal error: logical ind.list must have length nInd(x)\n")
    }
    idx <- which(ind.list)
  } else if (is.numeric(ind.list)) {
    if (any(ind.list < 1 | ind.list > n_ind)) {
      stop("Fatal error: numeric ind.list contains indices outside [1, nInd(x)]\n")
    }
    idx <- ind.list
  } else if (is.character(ind.list)) {
    idx <- match(ind.list, all_inds)
    if (any(is.na(idx))) {
      missing_inds <- ind.list[is.na(idx)]
      stop(
        "Fatal error: the following individuals in ind.list were not found in indNames(x):\n  ",
        paste(missing_inds, collapse = ", "),
        "\n"
      )
    }
  } else {
    stop("Fatal error: ind.list must be logical, numeric, or character\n")
  }
  
  if (length(idx) == 0) {
    warning("No individuals selected in ind.list; no changes made\n")
    return(x)
  }
  
  # ---- DO THE JOB ----
  
  current_pops <- pop(x)
  current_pops_char <- as.character(current_pops)
  
  old_pops_subset <- unique(current_pops_char[idx])
  
  # Assign new population label to selected individuals
  current_pops_char[idx] <- new.pop
  
  # Re-factor to ensure new.pop is included as a level
  pop(x) <- factor(current_pops_char)
  
  if (verbose >= 2) {
    cat(report(
      "  Reassigned",
      length(idx),
      "individual(s) from populations",
      paste(old_pops_subset, collapse = ", "),
      "to population",
      new.pop,
      "\n"
    ))
  }
  
  if (verbose >= 3) {
    cat("  Summary of recoded dataset\n")
    cat(paste("    No. of loci:       ", nLoc(x), "\n"))
    cat(paste("    No. of individuals:", nInd(x), "\n"))
    cat(paste("    No. of populations:", nPop(x), "\n"))
  }
  
  # ---- ADD TO HISTORY ----
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(x)
}
