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
.fbm_col_blocks <- function(p, nr) {
  size <- max(1L, as.integer(bigstatsr::block_size(max(1L, nr)) %/% 16L))
  starts <- seq(1L, p, by = size)
  lapply(starts, function(s) s:min(p, s + size - 1L))
}

## Applies FUN, a column-wise statistic of the decoded genotype matrix (one
## value per locus, e.g. colMeans(t == 1, na.rm = TRUE)), and returns the
## same vector as FUN(as.matrix(x)). For an FBM-backed object, FUN gets
## blocks of loci decoded exactly as by as.matrix() (integer storage for
## integer codes, indNames/locNames as dimnames), and the results are
## concatenated, so the whole matrix is never decoded and the values are
## identical. For a gen-backed object, FUN(as.matrix(x)).
.gl_colwise <- function(x, FUN) {
  G <- .fbm_or_null(x)
  if (is.null(G)) return(FUN(as.matrix(x)))
  code <- G$code256
  int <- all(code[!is.na(code)] %% 1 == 0)
  ind <- indNames(x)
  loc <- locNames(x)
  no_dn <- is.null(ind) && is.null(loc)
  out <- lapply(.fbm_col_blocks(ncol(G), nrow(G)), function(cols) {
    B <- G[, cols, drop = FALSE]
    if (int) storage.mode(B) <- "integer"
    if (!no_dn) dimnames(B) <- list(ind, loc[cols])
    FUN(B)
  })
  do.call(c, out)
}
