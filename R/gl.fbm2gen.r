#' Convert an FBM-backed dartR to a GEN-backed dartR (streamed with big_apply)
#'
#' @description
#' `gl.fbm2gen()` converts a `dartR` whose genotypes live in the `@fbm` slot into
#' a `dartR` with genotypes in the `@gen` (list of SNPbin) slot. The operation is
#' **column-chunked** via `bigstatsr::big_apply`: each block is decoded from FBM,
#' turned into a small `genlight`, and concatenated using `cbind.dartR` (SNPbin
#' path). At the end, `@fbm` is set to `NULL` and `@gen` holds the SNPbin list.
#'
#' @param x A `dartR` object with `@fbm` populated. If `@fbm` is `NULL`,
#'   the object is returned unchanged.
#' @param chunk Integer, number of **loci per block** to read from FBM
#'   (default `2048L`). Increase for speed if you have RAM to spare.
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
#' \dontrun{
#' d_gen <- gl.fbm2gen(d_fbm, chunk = 4096L)
#' length(d_gen@gen)      # nInd
#' nLoc(d_gen)            # number of loci (via genlight path)
#' }
#' @export

gl.fbm2gen <- function(x, chunk = 2048L, quiet = TRUE, verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  stopifnot(inherits(x, "dartR"))
  ## Safe FBM accessor (tolerates missing slot)
  .fbm_or_null <- function(obj) tryCatch(methods::slot(obj, "fbm"), error = function(e) NULL)

  fbm <- .fbm_or_null(x)
  if (is.null(fbm)) {
    if (!quiet && verbose >= 2) {
      message("gl.fbm2gen: no FBM found; returning input unchanged.")
    }
    return(x)
  }

  dummy <- new("genlight", gen=x@fbm[])
    
  x@gen <- dummy@gen
  x@fbm <- NULL
  ## Update locus-wise metadata that may have been concatenated during blocks
  
  methods::validObject(x)
  x
}

