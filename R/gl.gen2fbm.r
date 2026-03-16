#' Convert a gen-backed object to FBM-backed (streamed with big_apply)
#'
#' @description
#' `gl.gen2fbm()` converts a `dartR`/`genlight` object that stores genotypes
#' in the standard list-of-`SNPbin` format (`@gen`) into an `FBM.code256`
#' stored in the `@fbm` slot (and clears `@gen`, enforcing XOR).
#' The copy is done **by column blocks** via `bigstatsr::big_apply`, so the
#' full genotype matrix is never materialized in memory.
#'
#' @param x A **`dartR`** (preferred) or **`genlight`** object that currently
#'   has genotype data only in `@gen` (i.e., no FBM yet).
#' @param code A `bigsnpr` code mapping for `FBM.code256`. Defaults to
#'   `bigsnpr::CODE_012` (0/1/2 with `NA` support). 
#' @param backingfile File stem for the FBM backing files. Defaults to a temp file.
#' @param chunk Integer, number of **loci per block** to write with `big_apply`
#'   (default `2048L`). Increase for faster IO if you have RAM to spare.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' @return A **`dartR`** object with `@fbm` populated and `@gen` emptied.
#' @examples
#' \dontrun{
#' library(adegenet); library(bigstatsr); library(bigsnpr)
#' # say 'gl' is a genlight or dartR without FBM
#' d_fbm <- gl.gen2fbm(gl, code = bigsnpr::CODE_012, chunk = 4096L)
#' nInd(d_fbm); nLoc(d_fbm)   # dimensions via FBM
#' head(as.matrix(d_fbm))     # decodes from FBM
#' }
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
  
  if (datatype!="SNP") {
    stop("Only SNP data supported at this time.")
  }
  
  ## returns NULL if the 'fbm' slot is missing OR is NULL
 # .fbm_or_null <- function(x) {
  #  if (methods::.hasSlot(x, "fbm")) {
  #    val <- methods::slot(x, "fbm")
  #    return(if (is.null(val)) NULL else val)
  #  }
  #  NULL
#  }
  fbm <- .fbm_or_null(x)
  if (!is.null(fbm) && length(x@gen) == 0L) {
    if (verbose>2) message("Object already FBM-backed; returning as-is.")
    return(x)
  }
  
  ## Accept genlight too; coerce to dartR if needed (gen-only mode)
  if (is.null(fbm)) {
  
  ## Must have genotypes in @gen to convert
  if (length(x@gen) == 0L)
    stop("No genotypes found in @gen; nothing to convert.")
  
  n <- nInd(x); p <- nLoc(x)
  if (n == 0L || p == 0L)
    stop("Empty object (0 individuals or 0 loci).")
  
  ## Create destination FBM
  G <- bigstatsr::FBM.code256(n, p, code = code, backingfile = backingfile)
  
  ## Column-chunked write: big_apply over 'G' columns in blocks of 'chunk'
  write_block <- function(Y, ind, ind.col) {
    Xblk <- as.matrix(x[, ind.col, drop = FALSE])     # (n x |ind.col|)
    Xblk[is.na(Xblk)] <- 3                           # NA code]
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
  x
} else {
    stop("Object already FBM-backed; cannot convert.")
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(x)
}
