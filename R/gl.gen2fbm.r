#' @name gl.gen2fbm
#' @title Converts a gen-backed genlight object to a file-backed (FBM) one
#' @family data manipulation
#'
#' @description
#' `gl.gen2fbm()` converts a `dartR`/`genlight` object that stores genotypes
#' in the standard list-of-`SNPbin` format (`@gen`) into an `FBM.code256`
#' stored in the `@fbm` slot (and clears `@gen`, enforcing XOR).
#' The copy is done **by column blocks** via `bigstatsr::big_apply`, so the
#' full genotype matrix is never materialized in memory.
#'
#' @param x A **`dartR`** (preferred) or **`genlight`** object that currently
#'   has genotype data only in `@gen` (i.e., no FBM yet). A plain `genlight`
#'   is converted to `dartR` first. Only SNP data are supported.
#' @param code A `bigsnpr` code mapping for `FBM.code256`. Defaults to
#'   `bigsnpr::CODE_012` (0/1/2 with `NA` support).
#' @param backingfile File stem for the FBM backing files. Defaults to a temp file.
#' @param chunk Integer, number of **loci per block** to write with `big_apply`
#'   (default `2048L`). Increase for faster IO if you have RAM to spare.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#' @return A **`dartR`** object with `@fbm` populated and `@gen` emptied. An
#' object that is already FBM-backed is returned unchanged.
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' x <- gl.gen2fbm(testset.gl)
#' nInd(x); nLoc(x)            # dimensions via FBM
#' as.matrix(x)[1:5, 1:5]      # decodes from FBM
#' @seealso \code{\link{gl.fbm2gen}}
#' @export
gl.gen2fbm <- function(x,
                       code        = bigsnpr::CODE_012,
                       backingfile = tempfile("geno_"),
                       chunk       = 2048L,
                       verbose     = NULL) {

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2025.1",
                   verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (datatype != "SNP") {
    stop(error("  Only SNP data can be converted to FBM at this time\n"))
  }

  # A plain genlight has no @fbm slot; as() gives a valid dartR object
  # (class<- does not add the slot)
  if (!is(x, "dartR")) {
    x <- methods::as(x, "dartR")
  }

  fbm <- .fbm_or_null(x)
  if (!is.null(fbm)) {
    if (length(x@gen) > 0L) {
      stop(error("  Invalid object: both @fbm and @gen hold genotypes\n"))
    }
    if (verbose >= 2) {
      cat(report("  Object already FBM-backed; returned unchanged\n"))
    }
    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }
    return(x)
  }

  ## Must have genotypes in @gen to convert
  if (length(x@gen) == 0L) {
    stop(error("  No genotypes found in @gen; nothing to convert\n"))
  }

  n <- nInd(x); p <- nLoc(x)
  if (n == 0L || p == 0L) {
    stop(error("  Empty object (0 individuals or 0 loci)\n"))
  }

  # DO THE JOB

  ## Create destination FBM
  G <- bigstatsr::FBM.code256(n, p, code = code, backingfile = backingfile)

  ## Column-chunked write: big_apply over 'G' columns in blocks of 'chunk'
  write_block <- function(Y, ind, ind.col) {
    Xblk <- as.matrix(x[, ind.col, drop = FALSE])     # (n x |ind.col|)
    Xblk[is.na(Xblk)] <- 3                           # NA code
    Y[ind, ind.col] <- Xblk[ind, , drop = FALSE]
    NULL
  }

  bigstatsr::big_apply(
    G,
    a.FUN      = write_block,
    a.combine  = "c",
    ind        = seq_len(n),                          # explicit rows
    ind.col    = seq_len(p),                          # explicit columns
    block.size = max(1L, as.integer(chunk))           # columns per block
  )

  ## Set FBM and clear heavy SNPbin list (XOR)
  x@fbm <- G
  x@gen <- vector("list", 0L)

  methods::validObject(x)

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(x)
}
