#' @name gl.fbm.cleanup
#' @title Removes backing files of file-backed objects that are no longer used
#' @family data manipulation
#'
#' @description
#' File-backed (FBM) dartR objects (see \code{\link{gl.gen2fbm}}) keep their
#' genotypes in a backing file (.bk). Every subset, filter or imputation of a
#' file-backed object writes a new backing file. dartR deletes a backing file
#' automatically once no object in the R session uses it any more, when R's
#' garbage collector removes the last such object.
#'
#' \code{gl.fbm.cleanup()} forces that step: it runs the garbage collector,
#' deletes the backing files that are no longer used, and reports what is
#' left. Use it in long sessions, or before checking the free disk space.
#'
#' @details
#' Only backing files created by dartR in the current session are ever
#' deleted. A backing file named by the user (argument \code{backingfile} of
#' \code{\link{gl.gen2fbm}}, \code{cbind} or \code{rbind}) is never deleted.
#'
#' On Windows a backing file cannot be deleted while it is still
#' memory-mapped, which it can be for a short while after its last object is
#' removed. Such files are deleted at the next attempt (by this function, or
#' when dartR creates the next backing file).
#'
#' Backing files are placed in \code{getOption("FBM.dir")}, by default the
#' session temporary folder, which R deletes at the end of the session. If
#' \code{FBM.dir} is set to another folder, files that could not be deleted
#' when R closed can remain there; they can be deleted by hand once no R
#' session uses them.
#'
#' A file stays while anything in the session still refers to its object.
#' ggplot2 keeps the last plot drawn (for \code{ggsave()} and
#' \code{last_plot()}), and that plot refers to the input object of the dartR
#' function that drew it, so one earlier backing file can stay in use until
#' the next plot is drawn.
#'
#' Objects saved with \code{saveRDS()} and loaded again with
#' \code{readRDS()} in the same session are not known to dartR: their backing
#' file can be deleted when the original object is removed. Use
#' \code{\link{gl.save}} and \code{\link{gl.load}} instead.
#'
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#' @return Invisibly, a list with the number and size (Mb) of the backing
#' files deleted (\code{deleted}, \code{deleted.mb}), still in use
#' (\code{in.use}, \code{in.use.mb}) and still waiting to be deleted
#' (\code{pending}).
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' x <- gl.gen2fbm(testset.gl)
#' for (i in 1:5) x <- x[-1, ]   # each subset writes a new backing file
#' gl.fbm.cleanup()
#' @seealso \code{\link{gl.gen2fbm}}, \code{\link{gl.save}}
#' @export

gl.fbm.cleanup <- function(verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  funname <- utils.flag.start(func = funname, verbose = verbose)

  # DO THE JOB
  mb <- function(paths) {
    if (!length(paths)) return(0)
    sum(file.size(paths), na.rm = TRUE) / 1024^2
  }

  # Files that will be deleted below: pending ones, and those whose objects
  # the garbage collector is about to remove (counted after the run)
  before <- c(ls(.dartR_fbm_files, all.names = TRUE),
              ls(.dartR_fbm_pending, all.names = TRUE))
  before <- unique(before[file.exists(before)])
  size_before <- vapply(before, function(f) file.size(f), numeric(1))

  # Two runs: the first runs the finalizers of unused objects (which delete
  # their files or mark them pending), the second releases memory mappings
  invisible(gc(verbose = FALSE))
  invisible(gc(verbose = FALSE))
  .fbm_retry_pending()

  deleted <- before[!file.exists(before)]
  in_use <- ls(.dartR_fbm_files, all.names = TRUE)
  in_use <- in_use[file.exists(in_use)]
  pending <- ls(.dartR_fbm_pending, all.names = TRUE)

  res <- list(deleted = length(deleted),
              deleted.mb = sum(size_before[deleted]) / 1024^2,
              in.use = length(in_use),
              in.use.mb = mb(in_use),
              pending = length(pending))

  if (verbose >= 2) {
    cat(report("  Deleted", res$deleted, "unused backing file(s),",
               round(res$deleted.mb, 1), "Mb\n"))
    cat(report("  In use:", res$in.use, "backing file(s),",
               round(res$in.use.mb, 1), "Mb, in",
               getOption("FBM.dir", tempdir()), "\n"))
    if (res$pending > 0) {
      cat(warn("  Warning:", res$pending, "backing file(s) still locked by",
               "the system; they will be deleted at the next attempt\n"))
    }
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  invisible(res)
}
