## Creation of file-backed matrices (FBM) for dartR objects
##
## Every temporary backing file that dartR creates goes through these helpers,
## so that
##   - the file is placed in getOption("FBM.dir"), the folder bigstatsr uses
##     for its own backing files (by default the session temporary folder),
##   - the file is recorded as created by dartR. Only recorded files may ever
##     be removed automatically; a backing file named by the user (argument
##     backingfile = ...) is never recorded.

## Backing files created by dartR in this session (names: full .bk paths)
.dartR_fbm_files <- new.env(parent = emptyenv())

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

## Record a backing file as created by dartR (internal)
.fbm_register <- function(G) {
  assign(normalizePath(G$backingfile, winslash = "/", mustWork = FALSE),
         TRUE, envir = .dartR_fbm_files)
  invisible(G)
}

## TRUE if the backing file of G was created by dartR in this session
.fbm_is_registered <- function(G) {
  exists(normalizePath(G$backingfile, winslash = "/", mustWork = FALSE),
         envir = .dartR_fbm_files, inherits = FALSE)
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
  if (temp) backingfile <- utils.fbm.tempfile(prefix)
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
  G2 <- bigstatsr::big_copy(G, ind.row = ind.row, ind.col = ind.col,
                            backingfile = utils.fbm.tempfile(prefix))
  .fbm_register(G2)
}
