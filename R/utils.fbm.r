## Creation of file-backed matrices (FBM) for dartR objects
##
## Every temporary backing file that dartR creates goes through these helpers,
## so that
##   - the file is placed in getOption("FBM.dir"), the folder bigstatsr uses
##     for its own backing files (by default the session temporary folder),
##   - the file is recorded as created by dartR. Only recorded files may ever
##     be removed automatically; a backing file named by the user (argument
##     backingfile = ...) is never recorded.
##
## Automatic clean-up. Each R object (FBM) that points to a recorded file is
## counted, and a finalizer on the object lowers the count when the object
## is garbage collected. When the count reaches zero no object uses the file
## any more and it is deleted. Counting objects, not files, matters because
## some bigstatsr/bigsnpr functions (e.g. snp_fastImputeSimple) return a new
## R object on the same backing file. Copies of a dartR object share one FBM
## object and are covered by the same count.
##
## On Windows a file cannot be deleted while it is still memory-mapped, which
## it can be for a short while after the object is collected. Such files are
## kept in a pending list and deleted at the next attempt: when dartR creates
## the next backing file, in gl.fbm.cleanup(), and when the package unloads.

## Backing files created by dartR in this session: full .bk path -> number of
## live R objects that point to the file
.dartR_fbm_files <- new.env(parent = emptyenv())
## Files whose deletion failed (still mapped); retried later
.dartR_fbm_pending <- new.env(parent = emptyenv())

#' @name utils.fbm.tempfile
#' @title File stem for a new temporary FBM backing file
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' Returns a file stem (without the .bk extension) in the folder given by
#' \code{getOption("FBM.dir")}, the option bigstatsr uses for its own backing
#' files, or in \code{tempdir()} if the option is not set. The folder is
#' created if it does not exist. To keep file-backed genotypes on another
#' disk, set e.g. \code{options(FBM.dir = "D:/scratch")}.
#' @param prefix Start of the file name [default "geno_"].
#' @return A file stem.
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.fbm.tempfile <- function(prefix = "geno_") {
  dir <- getOption("FBM.dir")
  if (is.null(dir) || !nzchar(dir)) dir <- tempdir()
  if (!dir.exists(dir) && !dir.create(dir, recursive = TRUE)) {
    stop(error("  Fatal Error: cannot create the folder for FBM backing files",
               "(option FBM.dir):", dir, "\n"))
  }
  tempfile(prefix, tmpdir = dir)
}

.fbm_path <- function(G) {
  normalizePath(G$backingfile, winslash = "/", mustWork = FALSE)
}

## Delete a backing file (and a .rds descriptor, if any); if it is still
## mapped, keep it in the pending list. Never raises an error: it can run
## inside a finalizer.
.fbm_try_delete <- function(path) {
  tryCatch({
    unlink(c(path, sub("\\.bk$", ".rds", path)))
    if (file.exists(path)) {
      assign(path, TRUE, envir = .dartR_fbm_pending)
    } else if (exists(path, envir = .dartR_fbm_pending, inherits = FALSE)) {
      rm(list = path, envir = .dartR_fbm_pending)
    }
  }, error = function(e) NULL)
  invisible(NULL)
}

## Retry the deletion of pending files
.fbm_retry_pending <- function() {
  for (path in ls(.dartR_fbm_pending, all.names = TRUE)) .fbm_try_delete(path)
  invisible(NULL)
}

## Finalizer factory: the closure holds only the path, never the object
.fbm_release_fun <- function(path) {
  force(path)
  function(e) {
    tryCatch({
      cnt <- get0(path, envir = .dartR_fbm_files, inherits = FALSE,
                  ifnotfound = 0L) - 1L
      if (cnt > 0L) {
        assign(path, cnt, envir = .dartR_fbm_files)
      } else {
        if (exists(path, envir = .dartR_fbm_files, inherits = FALSE)) {
          rm(list = path, envir = .dartR_fbm_files)
        }
        .fbm_try_delete(path)
      }
    }, error = function(e) NULL)
  }
}

## Record an R object pointing to a backing file created by dartR, so that
## the file is deleted once no recorded object uses it (internal). Use it
## also for a new R object that bigstatsr/bigsnpr returns on a recorded file.
## Recording the same object twice has no effect.
.fbm_register <- function(G) {
  env <- as.environment(G)
  if (isTRUE(get0(".dartR_registered", envir = env, inherits = FALSE))) {
    return(invisible(G))
  }
  path <- .fbm_path(G)
  cnt <- get0(path, envir = .dartR_fbm_files, inherits = FALSE,
              ifnotfound = 0L)
  assign(path, cnt + 1L, envir = .dartR_fbm_files)
  assign(".dartR_registered", TRUE, envir = env)
  reg.finalizer(env, .fbm_release_fun(path), onexit = TRUE)
  invisible(G)
}

## TRUE if the backing file of G was created by dartR in this session and is
## still in use
.fbm_is_registered <- function(G) {
  exists(.fbm_path(G), envir = .dartR_fbm_files, inherits = FALSE)
}

#' @name utils.fbm.new
#' @title Creates an FBM.code256 for genotypes in a temporary backing file
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' Creates a new \code{bigstatsr::FBM.code256}. With \code{backingfile =
#' NULL} the file is placed by \code{\link{utils.fbm.tempfile}} and recorded
#' as created by dartR; a backing file named by the caller is used as given
#' and not recorded.
#' @param n Number of individuals (rows) [required].
#' @param p Number of loci (columns) [required].
#' @param code Code mapping of the FBM [default bigsnpr::CODE_012].
#' @param backingfile File stem of the backing file, or NULL for a temporary
#' file [default NULL].
#' @param prefix Start of the temporary file name [default "geno_"].
#' @return An FBM.code256.
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.fbm.new <- function(n, p,
                          code = bigsnpr::CODE_012,
                          backingfile = NULL,
                          prefix = "geno_") {
  temp <- is.null(backingfile)
  if (temp) {
    .fbm_retry_pending()
    backingfile <- utils.fbm.tempfile(prefix)
  }
  G <- bigstatsr::FBM.code256(n, p, code = code, backingfile = backingfile)
  if (temp) .fbm_register(G)
  G
}

#' @name utils.fbm.copy
#' @title Copies (part of) an FBM.code256 into a temporary backing file
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' A wrapper of \code{bigstatsr::big_copy()} whose backing file is placed by
#' \code{\link{utils.fbm.tempfile}} and recorded as created by dartR.
#' @param G An FBM.code256 [required].
#' @param ind.row Rows to copy [default all].
#' @param ind.col Columns to copy [default all].
#' @param prefix Start of the temporary file name [default "geno_"].
#' @return An FBM.code256 with the same code mapping.
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @export
utils.fbm.copy <- function(G,
                           ind.row = bigstatsr::rows_along(G),
                           ind.col = bigstatsr::cols_along(G),
                           prefix = "geno_") {
  .fbm_retry_pending()
  G2 <- bigstatsr::big_copy(G, ind.row = ind.row, ind.col = ind.col,
                            backingfile = utils.fbm.tempfile(prefix))
  .fbm_register(G2)
}

## Blocks of loci (columns) for reading 'nr' rows of an FBM. The block size
## is bigstatsr's (which follows options(bigstatsr.block.sizeGB), default
## 1 Gb of doubles) divided by 16, as a block is held in several copies
## (decoded values, logical masks, results): about 64 Mb per copy with the
## default. Reading blocks of loci follows the storage order of the backing
## file (column by column).
.fbm_block_size <- function(nr) {
  max(1L, as.integer(bigstatsr::block_size(max(1L, nr)) %/% 16L))
}

.fbm_col_blocks <- function(p, nr) {
  size <- .fbm_block_size(nr)
  starts <- seq(1L, p, by = size)
  lapply(starts, function(s) s:min(p, s + size - 1L))
}

## Applies FUN, a column-wise statistic of the decoded genotype matrix (one
## value per locus, e.g. colMeans(t == 1, na.rm = TRUE)), and returns the
## same vector as FUN(as.matrix(x)). For an FBM-backed object, FUN gets
## blocks of loci decoded exactly as by as.matrix() (integer storage for
## integer codes, indNames/locNames as dimnames), and the results are
## concatenated, so the whole matrix is never decoded and the values are
## identical. For a gen-backed object, FUN(as.matrix(x)). With 'rows' (row
## indices, e.g. the individuals of one population), FUN gets only those
## rows, as FUN(as.matrix(x)[rows, , drop = FALSE]), without copying the
## object. FUN may return a matrix with one column per locus (several
## statistics); the blocks are then bound with cbind().
.gl_colwise <- function(x, FUN, rows = NULL) {
  out <- .gl_colblocks(x, FUN, rows = rows)
  # several statistics per locus (a matrix with loci in columns) are bound
  # column-wise, single statistics concatenated
  if (is.matrix(out[[1]])) do.call(cbind, out) else do.call(c, out)
}

## The building block of .gl_colwise: applies FUN to blocks of loci of an
## FBM-backed object (decoded as by as.matrix(), optionally only 'rows') and
## returns the list of per-block results, for callers that combine them
## themselves (e.g. FUN returns a list of matrices, or counts to be added).
## For a gen-backed object, a list with the single result
## FUN(as.matrix(x)[rows, ]).
.gl_colblocks <- function(x, FUN, rows = NULL) {
  G <- .fbm_or_null(x)
  if (is.null(G)) {
    m <- as.matrix(x)
    return(list(FUN(if (is.null(rows)) m else m[rows, , drop = FALSE])))
  }
  code <- G$code256
  int <- all(code[!is.na(code)] %% 1 == 0)
  ind <- indNames(x)
  loc <- locNames(x)
  no_dn <- is.null(ind) && is.null(loc)
  if (!is.null(rows) && !no_dn && !is.null(ind)) ind <- ind[rows]
  lapply(.fbm_col_blocks(ncol(G), if (is.null(rows)) nrow(G) else
           length(rows)), function(cols) {
    B <- if (is.null(rows)) G[, cols, drop = FALSE] else
      G[rows, cols, drop = FALSE]
    if (int) storage.mode(B) <- "integer"
    if (!no_dn) dimnames(B) <- list(ind, loc[cols])
    FUN(B)
  })
}

## Applies FUN, a per-individual statistic that is a sum over loci (e.g.
## rowSums(is.na(t))), and returns the same vector as FUN(as.matrix(x)). For
## an FBM-backed object, FUN gets blocks of loci decoded as by as.matrix()
## and the per-block results are added; sums of counts are exact in double
## precision, so the values are identical. Only for statistics that are
## sums over loci. For a gen-backed object, FUN(as.matrix(x)).
.gl_rowwise <- function(x, FUN) {
  G <- .fbm_or_null(x)
  if (is.null(G)) return(FUN(as.matrix(x)))
  code <- G$code256
  int <- all(code[!is.na(code)] %% 1 == 0)
  ind <- indNames(x)
  loc <- locNames(x)
  no_dn <- is.null(ind) && is.null(loc)
  res <- NULL
  for (cols in .fbm_col_blocks(ncol(G), nrow(G))) {
    B <- G[, cols, drop = FALSE]
    if (int) storage.mode(B) <- "integer"
    if (!no_dn) dimnames(B) <- list(ind, loc[cols])
    r <- FUN(B)
    res <- if (is.null(res)) r else res + r
  }
  res
}

## Number of missing genotypes, as sum(is.na(as.matrix(x))), in blocks of
## loci for FBM-backed objects
.gl_n_na <- function(x) {
  sum(.gl_colwise(x, function(t) colSums(is.na(t))))
}

## Genotype counts per population and locus: for each population (in the
## order of seppop(), the levels of pop(x)) a 3 x nLoc matrix with the
## numbers of genotypes 0, 1 and 2 (colSums(m == k, na.rm = TRUE), double,
## named by locus), read from each population's rows in one pass over
## blocks of loci, so the populations are not copied and a file-backed
## object is not decoded as a whole.
.gl_pop_genotype_counts <- function(x) {
  lv <- levels(pop(x))
  rows <- lapply(lv, function(p) which(pop(x) == p))
  blocks <- .gl_colblocks(x, function(B) {
    lapply(rows, function(r) {
      m <- B[r, , drop = FALSE]
      rbind(colSums(m == 0, na.rm = TRUE), colSums(m == 1, na.rm = TRUE),
            colSums(m == 2, na.rm = TRUE))
    })
  })
  res <- lapply(seq_along(lv), function(g) {
    do.call(cbind, lapply(blocks, function(b) b[[g]]))
  })
  names(res) <- lv
  res
}

## R's mean() of a logical or integer vector, rowMeans() and colMeans()
## divide a long-double sum by the count and round the quotient to double.
## For a statistic counted over blocks of loci, a plain num / den differs
## from that by one unit in the last place for about one count in 6,000
## once den exceeds 2,048. Returns, element by element, exactly the value
## that mean() (rowMeans(), colMeans()) gives for num TRUE among den values,
## so counted statistics stay identical to the decoded-matrix versions.
## num and den are whole numbers; den = 0 gives NaN, as rowMeans(na.rm =
## TRUE) does for a row without values.
.gl_count_mean <- function(num, den) {
  num <- as.numeric(num)
  den <- rep_len(as.numeric(den), length(num))
  res <- num / den
  ok <- which(!is.na(num) & !is.na(den) & den > 2048 & num > 0 & num < den)
  if (length(ok) > 0) {
    key <- paste(num[ok], den[ok])
    u <- !duplicated(key)
    val <- mapply(function(a, b) mean(rep(c(TRUE, FALSE), c(a, b - a))),
                  num[ok][u], den[ok][u])
    res[ok] <- val[match(key, key[u])]
  }
  names(res) <- names(num)
  res
}

## Pairwise genotype comparisons of all individuals of an FBM-backed object,
## counted over blocks of loci read from the backing file: for each pair of
## individuals (i, j), the number of loci scored in both (L), and among those
## the numbers that differ by one (n1: homozygote vs heterozygote) and by two
## allele copies (n2: opposite homozygotes). Assumes genotypes coded 0, 1,
## 2 (number of alternative alleles, diploid) and NA for missing, as dartR
## writes them (CODE_012). The counts are integers, so they do not depend
## on the block size. Returns a list of three integer vectors in the order
## of a 'dist' object (lower triangle by columns), each of length
## n(n - 1)/2: the pairwise counts take 12 bytes per pair, e.g. 2.4 Gb for
## 20,000 individuals, besides the n x n result the callers build.
.gl_pair_counts <- function(x) {
  G <- .fbm_or_null(x)
  if (is.null(G)) {
    stop(error("  Fatal Error: .gl_pair_counts() needs a file-backed object\n"))
  }
  code <- G$code256
  if (!all(code[!is.na(code)] %in% 0:2)) {
    stop(error("  Fatal Error: pairwise distances from a file-backed object",
               "need genotypes coded 0, 1, 2 or NA\n"))
  }
  # Stub to satisfy package checks; the real function is compiled below
  dartR_pair_counts_acc <- function(...) NULL
  Rcpp::cppFunction(
'
    void dartR_pair_counts_acc(const NumericMatrix& B, IntegerVector L,
                               IntegerVector n1, IntegerVector n2) {
      // Adds the comparisons of one block of loci (individuals x loci,
      // coded 0, 1, 2 or NA) to the pair counts L, n1, n2 (dist order),
      // in place.
      const int n = B.nrow(), b = B.ncol();
      // individuals in rows of a transposed copy, so that the loci of one
      // individual are contiguous; -1 marks a missing genotype
      std::vector<int> g((size_t) n * b);
      for (int k = 0; k < b; ++k) {
        for (int i = 0; i < n; ++i) {
          double v = B(i, k);
          g[(size_t) i * b + k] = ISNAN(v) ? -1 : (int) v;
        }
      }
      R_xlen_t idx = 0;
      for (int j = 0; j < n - 1; ++j) {
        const int* gj = &g[(size_t) j * b];
        for (int i = j + 1; i < n; ++i, ++idx) {
          const int* gi = &g[(size_t) i * b];
          int l = 0, d1 = 0, d2 = 0;
          for (int k = 0; k < b; ++k) {
            int a = gi[k], c = gj[k];
            if (a < 0 || c < 0) continue;
            ++l;
            int d = a > c ? a - c : c - a;
            d1 += (d == 1);
            d2 += (d == 2);
          }
          L[idx] += l;
          n1[idx] += d1;
          n2[idx] += d2;
        }
      }
    }'
  )
  n <- nrow(G)
  npair <- n * (n - 1) / 2
  # fresh vectors, filled in place by the compiled function
  L <- integer(npair)
  n1 <- integer(npair)
  n2 <- integer(npair)
  if (npair > 0) {
    for (cols in .fbm_col_blocks(ncol(G), n)) {
      dartR_pair_counts_acc(G[, cols, drop = FALSE], L, n1, n2)
    }
  }
  list(L = L, n1 = n1, n2 = n2)
}

## A symmetric n x n matrix from a vector in 'dist' order, with the given
## diagonal and no dimnames
.gl_dist_to_matrix <- function(v, n, diag = 0) {
  m <- matrix(diag, n, n)
  lt <- lower.tri(m)
  m[lt] <- v
  m <- t(m)
  m[lt] <- v
  m
}

## Reader of the decoded genotypes of x: returns function(rows = NULL,
## cols = NULL) giving as.matrix(x)[rows, cols, drop = FALSE] (NULL: all
## rows or all columns), with the dimnames and storage of as.matrix(). For
## an FBM-backed object only the requested part is read from the backing
## file, so callers that ask for blocks of loci (.fbm_col_blocks), single
## individuals or selected loci never decode the whole matrix. Reading one
## individual touches the whole backing file (stored by loci), so callers
## should read few rows. For a gen-backed object the matrix is decoded once
## and kept by the reader, as as.matrix(x) would.
.gl_reader <- function(x) {
  G <- .fbm_or_null(x)
  if (is.null(G)) {
    m <- as.matrix(x)
    return(function(rows = NULL, cols = NULL) {
      if (is.null(rows) && is.null(cols)) return(m)
      if (is.null(rows)) return(m[, cols, drop = FALSE])
      if (is.null(cols)) return(m[rows, , drop = FALSE])
      m[rows, cols, drop = FALSE]
    })
  }
  code <- G$code256
  int <- all(code[!is.na(code)] %% 1 == 0)
  ind <- indNames(x)
  loc <- locNames(x)
  no_dn <- is.null(ind) && is.null(loc)
  function(rows = NULL, cols = NULL) {
    rows <- if (is.null(rows)) bigstatsr::rows_along(G) else
      seq_len(nrow(G))[rows]
    cols <- if (is.null(cols)) bigstatsr::cols_along(G) else
      seq_len(ncol(G))[cols]
    B <- G[rows, cols, drop = FALSE]
    if (int) storage.mode(B) <- "integer"
    if (!no_dn) dimnames(B) <- list(ind[rows], loc[cols])
    B
  }
}

## PLINK binary fileset (.bed/.bim/.fam) read without decoding all genotypes
## in memory (gl.read.PLINK).
##
## .bim and .fam are read exactly as snpStats::read.plink reads them
## (read.table, na.strings c("0", "-9"), the same column names and row
## names), so the metadata are identical to the snpStats route.
.plink_read_meta <- function(stub, na.strings = c("0", "-9"), sep = ".") {
  df.bim <- utils::read.table(paste0(stub, ".bim"), comment.char = "",
                              as.is = TRUE, na.strings = na.strings)
  snps <- as.character(df.bim[, 2])
  if (ncol(df.bim) == 6) {
    names(df.bim) <- c("chromosome", "snp.name", "cM", "position",
                       "allele.1", "allele.2")
  } else {
    warning("non-standard .bim file")
  }
  rownames(df.bim) <- snps
  if (any(duplicated(snps))) stop("duplicated snp name(s)")
  df.fam <- utils::read.table(paste0(stub, ".fam"), comment.char = "",
                              as.is = TRUE, na.strings = na.strings)
  ped <- as.character(df.fam[, 1])
  mem <- as.character(df.fam[, 2])
  if (any(duplicated(ped))) {
    if (any(duplicated(mem))) {
      id <- paste(ped, mem, sep = sep)
      if (any(duplicated(id))) {
        stop("couldn't create unique subject identifiers")
      }
    } else {
      id <- mem
    }
  } else {
    id <- ped
  }
  names(df.fam) <- c("pedigree", "member", "father", "mother", "sex",
                     "affected")
  rownames(df.fam) <- id
  list(fam = df.fam, map = df.bim)
}

## Genotypes of a SNP-major .bed file written in blocks of SNPs into a new
## temporary FBM.code256 (CODE_012). In a .bed each SNP takes ceiling(n/4)
## bytes, four individuals per byte from the low bits up, with the 2-bit
## codes 00 = homozygous allele.1, 01 = missing, 10 = heterozygous,
## 11 = homozygous allele.2 (allele.1/allele.2 as columns 5/6 of the .bim).
## As in snpStats (and gl.read.PLINK before), genotypes are stored as the
## count of allele.2: 0, 1, 2, and missing as the raw code 3 (NA).
## Also returns, per individual, the numbers of called and heterozygous
## genotypes, counted in the same pass. Only one block is in memory.
.plink_bed_to_fbm <- function(bedfile, n, p) {
  con <- file(bedfile, open = "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, "raw", 3L)
  if (length(magic) < 3L || magic[1] != as.raw(0x6c) ||
      magic[2] != as.raw(0x1b)) {
    stop(error("  Fatal Error:", bedfile, "is not a PLINK .bed file\n"))
  }
  if (magic[3] != as.raw(0x01)) {
    stop(error("  Fatal Error:", bedfile, "is in individual-major mode;",
               "convert it to SNP-major with plink --make-bed\n"))
  }
  nb <- (n + 3L) %/% 4L
  expected <- 3 + as.numeric(nb) * p
  if (file.size(bedfile) != expected) {
    stop(error("  Fatal Error: size of", bedfile, "does not match the",
               n, "individuals and", p, "SNPs of the .fam and .bim files\n"))
  }
  G <- utils.fbm.new(n, p, code = bigsnpr::CODE_012, prefix = "plink_")
  # 2-bit code -> stored value: 00 -> 0, 01 -> 3 (NA), 10 -> 1, 11 -> 2
  lut <- c(0L, 3L, 1L, 2L)
  called <- het <- numeric(n)
  for (cols in .fbm_col_blocks(p, n)) {
    v <- as.integer(readBin(con, "raw", nb * length(cols)))
    # codes of the four individuals of each byte, in individual order
    codes <- rbind(v %% 4L, (v %/% 4L) %% 4L, (v %/% 16L) %% 4L,
                   (v %/% 64L) %% 4L)
    codes <- matrix(codes, nrow = 4L * nb)[seq_len(n), , drop = FALSE]
    called <- called + rowSums(codes != 1L)
    het <- het + rowSums(codes == 2L)
    G[, cols] <- lut[codes + 1L]
  }
  list(G = G, called = called, het = het)
}
