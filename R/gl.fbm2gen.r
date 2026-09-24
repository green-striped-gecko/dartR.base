#' @name gl.fbm2gen
#' @title Converts a file-backed (FBM) genlight object to a gen-backed one
#' @family data manipulation
#'
#' @description
#' `gl.fbm2gen()` converts a `dartR` whose genotypes live in the `@fbm` slot
#' into a `dartR` with genotypes in the `@gen` (list of SNPbin) slot. The
#' conversion runs in blocks of individuals: each block is decoded from the
#' FBM and turned into SNPbin objects, so the full genotype matrix is never
#' decoded at once. At the end, `@fbm` is set to `NULL` and `@gen` holds the
#' SNPbin list.
#'
#' @param x A `dartR` or `genlight` object. If it holds no FBM, it is returned
#'   unchanged.
#' @param chunk Integer, number of **individuals per block** decoded from the
#'   FBM at a time; a block holds `chunk` x nLoc(x) values in memory
#'   (default `256L`). Increase for speed if you have RAM to spare.
#' @param quiet Logical; if `TRUE`, suppress the non-critical
#'   "no FBM found" message regardless of verbosity [default TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'   progress log; 3, progress and results summary; 5, full report
#'   [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#'   or 2 if no global is set].
#'
#' @return A `dartR` object with **`@gen` populated** and **`@fbm = NULL`**.
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' x <- gl.gen2fbm(testset.gl)
#' y <- gl.fbm2gen(x)
#' length(y@gen)      # one SNPbin per individual
#' identical(as.matrix(y), as.matrix(testset.gl))
#' @seealso \code{\link{gl.gen2fbm}}
#' @export

gl.fbm2gen <- function(x, chunk = 256L, quiet = TRUE, verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # A plain genlight cannot hold an FBM: nothing to convert
  fbm <- .fbm_or_null(x)
  if (is.null(fbm)) {
    if (!quiet && verbose >= 2) {
      cat(report("  No FBM found; returning input unchanged\n"))
    }
    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }
    return(x)
  }

  # DO THE JOB
  # Decode blocks of individuals so that the full matrix is never held in
  # memory; each row of the FBM becomes one SNPbin
  n <- nrow(fbm)
  chunk <- max(1L, as.integer(chunk))
  gen <- vector("list", n)
  for (s in seq(1L, n, by = chunk)) {
    rows <- s:min(n, s + chunk - 1L)
    blk <- fbm[rows, , drop = FALSE]
    gen[rows] <- methods::new("genlight", gen = blk,
                              ploidy = x@ploidy[rows])@gen
  }

  x@gen <- gen
  x@fbm <- NULL

  methods::validObject(x)

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  x
}
