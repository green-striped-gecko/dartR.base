#' @name gl.subsample.loci
#' @title Subsamples n loci from a genlight object [deprecated]
#' @family data manipulation
#'
#' @description
#' \strong{Deprecated.} Use \code{gl.subsample.loc(x, n, replace = FALSE,
#' method = method, mono.rm = mono.rm)} instead; see
#' \code{\link{gl.subsample.loc}}.
#'
#' Subsamples a genlight object on loci, at random or on information
#' content, by calling \code{gl.subsample.loc}. It will be removed in a
#' future release.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#'  (SilicoDArT) data [required].
#' @param n Number of loci to include in the subsample [required].
#' @param method Method: 'random', in which case the loci are sampled at random;
#' or 'pic', in which case the top n loci ranked on information content are
#' chosen. Information content is AvgPIC in the case of SNP data and PIC in
#' the case of presence/absence (SilicoDArT) data, recalculated first if the
#' stored values are out of date [default 'random'].
#' @param mono.rm Delete monomorphic loci before sampling [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return A genlight object with n loci
#' @export
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @seealso \code{\link{gl.subsample.loc}}
#' @examples
#' \donttest{
#' # SNP data
#' gl2 <- gl.subsample.loci(testset.gl, n=200, method='pic')
#' # Tag P/A data
#' gl2 <- gl.subsample.loci(testset.gs, n=100, method='random')
#' }

gl.subsample.loci <- function(x,
                              n,
                              method = "random",
                              mono.rm = FALSE,
                              verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  .Deprecated(
    new = "gl.subsample.loc",
    msg = paste0(
      "gl.subsample.loci() is deprecated and will be removed in a future ",
      "release.\nUse gl.subsample.loc(x, n, replace = FALSE, method = '",
      tolower(method), "', mono.rm = ", mono.rm, ") instead."
    )
  )

  # Unlike gl.subsample.loc, which caps n, this function has always stopped
  # on an out-of-range n; keep that behaviour
  n.avail <- if (mono.rm) {
    nLoc(gl.filter.monomorphs(x, verbose = 0))
  } else {
    nLoc(x)
  }
  if (missing(n) || n <= 0 || n > n.avail) {
    stop(error("  Fatal Error: subsample size must be a positive integer",
               ">= 1 and <=", n.avail, "\n"))
  }

  x.new <- gl.subsample.loc(x,
                            n = n,
                            replace = FALSE,
                            method = method,
                            mono.rm = mono.rm,
                            verbose = 0)

  # Record this call rather than the inner gl.subsample.loc() call
  nh <- length(x.new@other$history)
  x.new@other$history[[nh]] <- match.call()

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(x.new)
}
