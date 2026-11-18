#' @name gl.fbm2gen
#' @title Converts a file-backed (FBM) genlight object to a gen-backed one
#' @family data manipulation
#'
#' @description
#' `gl.fbm2gen()` converts a `dartR` whose genotypes live in the `@fbm` slot
#' into a `dartR` with genotypes in the `@gen` (list of SNPbin) slot. At the
#' end, `@fbm` is set to `NULL` and `@gen` holds the SNPbin list.
#'
#' The conversion reads the FBM in **blocks of loci** (the order in which the
#' backing file is stored, so reading is fast) and packs each block directly
#' into the bits of the SNPbin objects; the full genotype matrix is never
#' decoded at once. The block size follows the bigstatsr option
#' \code{bigstatsr.block.sizeGB} (default 1): a block holds one sixteenth of
#' \code{bigstatsr::block_size(nInd(x))} loci (at least 512). Beyond the
#' result itself, the conversion needs a few hundred Mb with the default
#' (about 370 Mb for 2,000 individuals x 200,000 loci); lower the option to
#' convert with less memory.
#'
#' Genotypes must be 0, 1, 2 or missing (0 or 1 for haploid individuals);
#' other values give an error.
#'
#' @param x A `dartR` or `genlight` object. If it holds no FBM, it is returned
#'   unchanged.
#' @param chunk Deprecated and ignored: the blocks are now blocks of loci
#'   whose size follows \code{bigstatsr.block.sizeGB}. Kept for
#'   compatibility [default NULL].
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

gl.fbm2gen <- function(x, chunk = NULL, quiet = TRUE, verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  funname <- utils.flag.start(func = funname, verbose = verbose)

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
  if (!is.null(chunk) && verbose >= 1) {
    cat(warn("  Argument 'chunk' is deprecated and ignored; the block size",
             "is set by options(bigstatsr.block.sizeGB)\n"))
  }

  ## Loci per block, a multiple of 8 so that blocks start on a byte of the
  ## packed SNPbin bits. A block is read as doubles and held in several
  ## copies (read, transposed, logical masks), hence 1/16 of bigstatsr's
  ## block size; at least 512 loci, as many small blocks add overhead
  loc_per <- max(512L, 8L * (bigstatsr::block_size(nrow(fbm)) %/% 128L))

  gen <- utils.fbm2gen.copy(fbm, ploidy = x@ploidy, loc_per = loc_per)

  x@gen <- gen
  x@fbm <- NULL

  methods::validObject(x)

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  x
}

## Converts the FBM G (n x p, decoded values 0, 1, 2, NA) into a list of
## SNPbin objects (one per individual), in blocks of 'loc_per' loci (a
## multiple of 8).
##
## SNPbin layout (adegenet): @snp holds bit planes (raw vectors), with bit k
## of byte b (least significant bit first) for locus 8 * (b - 1) + k.
## Missing genotypes have no bit set and are listed in @NA.posi.
## The planes follow adegenet's constructor (new("SNPbin", <integer>), as
## used by new("genlight", <matrix>)), which sets them from the individual's
## observed genotypes, not from the genlight ploidy: an individual with at
## least one genotype 2 gets two planes (genotype == 2, genotype >= 1) and
## SNPbin ploidy 2; otherwise one plane (genotype == 1) and SNPbin ploidy 1;
## an individual with all genotypes missing gets one empty plane and SNPbin
## ploidy NA. The result is identical to the SNPbins that the previous
## gl.fbm2gen (and new("genlight", <matrix>)) created.
utils.fbm2gen.copy <- function(G, ploidy, loc_per) {
  n <- nrow(G)
  p <- ncol(G)
  ploidy <- as.integer(rep_len(ploidy, n))
  if (anyNA(ploidy) || any(ploidy < 1L | ploidy > 2L)) {
    stop(error("  Only haploid and diploid data can be converted from FBM\n"))
  }
  hap <- which(ploidy == 1L)
  nbytes <- ceiling(p / 8)
  bytes_per <- max(1L, as.integer(loc_per %/% 8L))
  starts <- seq(1L, nbytes, by = bytes_per)

  ge1 <- rep(list(raw(nbytes)), n)            # bit set where genotype >= 1
  ge2 <- rep(list(raw(nbytes)), n)            # bit set where genotype == 2
  gna <- rep(list(raw(nbytes)), n)            # bit set where missing
  has2 <- logical(n)                          # any genotype 2 observed
  bad <- function() {
    stop(error("  Genotypes other than 0, 1, 2 (or 0, 1 for haploid",
               "individuals) found; they cannot be stored as SNPbin\n"))
  }

  for (k in seq_along(starts)) {
    b1 <- starts[k]
    b2 <- min(nbytes, b1 + bytes_per - 1L)
    b <- b1:b2
    j1 <- (b1 - 1L) * 8L + 1L
    j2 <- min(p, b2 * 8L)

    g <- t(G[, j1:j2, drop = FALSE])           # loci x n, decoded
    pad <- 8L * length(b) - nrow(g)
    if (pad > 0L) g <- rbind(g, matrix(0, pad, n))
    isna <- is.na(g)
    pna <- packBits(isna, "raw")
    g[isna] <- 0
    rm(isna)
    if (any(g != 0 & g != 1 & g != 2)) bad()
    p1 <- packBits(g >= 1, "raw")
    p2 <- packBits(g == 2, "raw")
    rm(g)
    dim(p1) <- dim(p2) <- dim(pna) <- c(length(b), n)
    has2_blk <- colSums(p2 != as.raw(0L)) > 0
    if (length(hap) && any(has2_blk[hap])) bad()
    has2 <- has2 | has2_blk
    for (i in seq_len(n)) {
      ge1[[i]][b] <- p1[, i]
      ge2[[i]][b] <- p2[, i]
      gna[[i]][b] <- pna[, i]
    }
    rm(p1, p2, pna)
  }

  proto <- methods::new("SNPbin")
  gen <- vector("list", n)
  for (i in seq_len(n)) {
    e <- proto
    e@snp <- if (has2[i]) list(ge2[[i]], ge1[[i]]) else list(ge1[[i]])
    e@n.loc <- as.integer(p)
    e@NA.posi <- which(as.logical(rawToBits(gna[[i]])))
    e@ploidy <- if (length(e@NA.posi) == p) {
      NA_integer_
    } else if (has2[i]) 2L else 1L
    gen[[i]] <- e
    ge1[i] <- list(NULL)
    ge2[i] <- list(NULL)
    gna[i] <- list(NULL)
  }
  gen
}
