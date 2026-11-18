#' @name gl.gen2fbm
#' @title Converts a gen-backed genlight object to a file-backed (FBM) one
#' @family data manipulation
#'
#' @description
#' `gl.gen2fbm()` converts a `dartR`/`genlight` object that stores genotypes
#' in the standard list-of-`SNPbin` format (`@gen`) into an `FBM.code256`
#' stored in the `@fbm` slot (and clears `@gen`, enforcing XOR).
#'
#' The genotypes are copied in **blocks of loci** read directly from the
#' packed bits of the `SNPbin` objects, so the full genotype matrix is never
#' decoded in memory. Each value is held as one byte while it is copied. The
#' block size follows the bigstatsr option \code{bigstatsr.block.sizeGB}
#' (default 1): a block holds as many loci as
#' \code{bigstatsr::block_size(nInd(x))}, and the peak memory use of the
#' conversion is roughly that option in Gb. Lower it (e.g.
#' \code{options(bigstatsr.block.sizeGB = 0.1)}) to convert with less memory.
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
#' @param chunk Deprecated. Number of loci per block; if given, it replaces
#'   the block size derived from \code{bigstatsr.block.sizeGB} (rounded up to
#'   a multiple of 8). Kept for compatibility [default NULL].
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
                       chunk       = NULL,
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

  ## Loci per block, a multiple of 8 so that blocks start on a byte of the
  ## packed SNPbin bits
  if (!is.null(chunk)) {
    if (verbose >= 1) {
      cat(warn("  Argument 'chunk' is deprecated; the block size is set by",
               "options(bigstatsr.block.sizeGB)\n"))
    }
    loc_per <- max(1L, as.integer(chunk))
  } else {
    loc_per <- bigstatsr::block_size(n)
  }
  loc_per <- min(8L * ceiling(p / 8), 8L * ceiling(loc_per / 8))

  ## Create destination FBM
  G <- utils.fbm.new(n, p, code = code, backingfile = backingfile)

  ## On error remove the partly written backing file
  ok <- FALSE
  on.exit(if (!ok) unlink(c(G$backingfile, G$rds)), add = TRUE)

  utils.gen2fbm.copy(x@gen, G, loc_per = loc_per)
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


## Copies a list of SNPbin objects (one per individual) into the FBM G (n x p),
## in blocks of 'loc_per' loci (a multiple of 8).
##
## SNPbin layout (adegenet): @snp is a list with one raw vector per allele
## copy (length = ploidy); bit k of byte b (least significant bit first) is
## that copy's allele at locus 8 * (b - 1) + k. The genotype (count of the
## second allele) is the sum over copies; @NA.posi holds the loci with a
## missing genotype. With at most two copies a and b the sum is
## xor(a, b) + 2 * (a & b), computed on raw bytes (one byte per value).
## Missing genotypes are written as raw 3, the NA code of FBM.code256.
## The first individual is checked against adegenet's own decoding, so that a
## change of the SNPbin layout gives an error instead of wrong genotypes.
utils.gen2fbm.copy <- function(gen, G, loc_per) {
  n <- length(gen)
  p <- ncol(G)
  nbytes <- ceiling(p / 8)
  bytes_per <- max(1L, as.integer(loc_per %/% 8L))

  plane <- function(e, k, b) {
    if (length(e@snp) >= k) e@snp[[k]][b] else raw(length(b))
  }
  nap <- lapply(gen, function(e) e@NA.posi)
  if (any(vapply(nap, is.unsorted, logical(1)))) {
    nap <- lapply(nap, sort)
  }
  starts <- seq(1L, nbytes, by = bytes_per)
  ends <- pmin(p, (starts + bytes_per - 1L) * 8L)       # last locus per block
  ## position in NA.posi where each block ends, per individual (computed
  ## once: findInterval copies its table to double on each call)
  cuts <- vapply(nap, function(na) findInterval(c(0.5, ends + 0.5), na),
                 integer(length(ends) + 1L))
  if (!is.matrix(cuts)) cuts <- matrix(cuts, ncol = n)
  na3 <- as.raw(3L)
  row1 <- as.integer(gen[[1]])
  row1[is.na(row1)] <- 3L

  for (k in seq_along(starts)) {
    b1 <- starts[k]
    b2 <- min(nbytes, b1 + bytes_per - 1L)
    b <- b1:b2
    j1 <- (b1 - 1L) * 8L + 1L
    j2 <- min(p, b2 * 8L)
    keep <- seq_len(j2 - j1 + 1L)

    ## bits of each allele copy: (8 * bytes) x n, one value per locus
    a1 <- rawToBits(vapply(gen, plane, raw(length(b)), k = 1L, b = b))
    a2 <- rawToBits(vapply(gen, plane, raw(length(b)), k = 2L, b = b))
    blk <- xor(a1, a2) | rawShift(a1 & a2, 1L)
    rm(a1, a2)
    dim(blk) <- c(8L * length(b), n)
    blk <- t(blk[keep, , drop = FALSE])                 # n x loci, raw

    ## missing genotypes in this block (NA.posi is sorted)
    for (i in which(cuts[k + 1L, ] > cuts[k, ])) {
      blk[i, nap[[i]][(cuts[k, i] + 1L):cuts[k + 1L, i]] - j1 + 1L] <- na3
    }

    if (b1 == 1L &&
        !identical(as.integer(blk[1, ]), row1[keep])) {
      stop(error("  Unexpected genotype storage in the genlight object",
                 "(SNPbin layout); please report this to the dartR team\n"))
    }

    G[, j1:j2] <- blk
    rm(blk)
  }
  invisible(G)
}
