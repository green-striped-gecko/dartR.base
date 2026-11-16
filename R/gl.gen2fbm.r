#' @name gl.gen2fbm
#' @title Converts a gen-backed genlight object to a file-backed (FBM) one
#' @family data manipulation
#'
#' @description
#' `gl.gen2fbm()` converts a `dartR`/`genlight` object that stores genotypes
#' in the standard list-of-`SNPbin` format (`@gen`) into an `FBM.code256`
#' stored in the `@fbm` slot (and clears `@gen`, enforcing XOR).
#' The copy is done **by column blocks** via `bigstatsr::big_apply`.
#'
#' Only diploid SNP data (genotypes 0, 1, 2 or missing) can be stored. The
#' FBM uses the value 3 for a missing genotype, so data with higher dosages
#' (e.g. polyploids) are refused with an error instead of being silently
#' turned into missing values.
#'
#' **Saving.** The FBM holds only a reference to its backing file, which by
#' default is a temporary file deleted at the end of the R session. Save a
#' file-backed object with \code{\link{gl.save}} (and reload it with
#' \code{\link{gl.load}}), never with \code{saveRDS()}: \code{saveRDS()} stores
#' only the reference, not the genotypes. After saving, check that the file
#' size makes sense for the number of individuals and loci (a file of a few
#' Kb for thousands of loci means the genotypes were not saved).
#'
#' @param x A **`dartR`** (preferred) or **`genlight`** object that currently
#'   has genotype data only in `@gen` (i.e., no FBM yet). A plain `genlight`
#'   is converted to `dartR` first. Only SNP data are supported.
#' @param code A `bigsnpr` code mapping for `FBM.code256`. Defaults to
#'   `bigsnpr::CODE_012` (0/1/2 with `NA` support).
#' @param backingfile File stem for the FBM backing file. If NULL, a temporary
#'   file in the folder \code{getOption("FBM.dir")} (see
#'   \code{\link{utils.fbm.tempfile}}). A backing file named here is never
#'   removed automatically
#'   [default NULL].
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
                       backingfile = NULL,
                       chunk       = 2048L,
                       verbose     = NULL) {

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  funname <- utils.flag.start(func = funname,
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

  ## Only diploid data fit the 0/1/2 + NA(=3) coding: a dosage of 3 would be
  ## stored as the NA code and silently lost
  if (any(adegenet::ploidy(x) > 2L, na.rm = TRUE)) {
    stop(error("  Only diploid data can be converted to FBM: ploidy above 2",
               "found. Dosages of 3 or more cannot be stored (3 is the code",
               "for a missing genotype)\n"))
  }

  # DO THE JOB

  ## Create destination FBM
  G <- utils.fbm.new(n, p, code = code, backingfile = backingfile)

  ## Column-chunked write: big_apply over 'G' columns in blocks of 'chunk'
  write_block <- function(Y, ind, ind.col) {
    Xblk <- as.matrix(x[, ind.col, drop = FALSE])     # (n x |ind.col|)
    if (any(Xblk > 2, na.rm = TRUE)) {
      stop(error("  Genotypes above 2 found (e.g. polyploid dosages); they",
                 "cannot be stored in an FBM, where 3 codes a missing",
                 "genotype\n"))
    }
    Xblk[is.na(Xblk)] <- 3                           # NA code
    Y[ind, ind.col] <- Xblk[ind, , drop = FALSE]
    NULL
  }

  ## On error remove the partly written backing file
  ok <- FALSE
  on.exit(if (!ok) unlink(c(G$backingfile, G$rds)), add = TRUE)

  bigstatsr::big_apply(
    G,
    a.FUN      = write_block,
    a.combine  = "c",
    ind        = seq_len(n),                          # explicit rows
    ind.col    = seq_len(p),                          # explicit columns
    block.size = max(1L, as.integer(chunk))           # columns per block
  )
  ok <- TRUE

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
