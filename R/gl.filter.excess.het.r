#' @name gl.filter.excess.het
#' @title Filters excessively-heterozygous loci from a genlight object [deprecated]
#' @family matched filter

#' @description
#' \strong{Deprecated.} Use
#' \code{gl.filter.hwe(direction = 'excess', test.type = 'ChiSquare',
#' mult.comp.adj = TRUE, mult.comp.adj.method = 'fdr')} instead — see
#' \code{\link{gl.filter.hwe}}.
#'
#' This wrapper reproduces the workflow of Robledo-Ruiz et al. (2023) for
#' removing loci with a significant excess of heterozygotes (candidate
#' sex-linked or paralogous loci) by calling gl.filter.hwe with the
#' equivalent settings, and will be removed in a future release.

#' @details
#' The published workflow is reproduced via gl.filter.hwe(direction =
#' 'excess', min.hobs = 0.5, ...): only loci with observed heterozygosity of
#' at least 0.5 are tested, and the false-discovery-rate adjustment is
#' confined to that screened set, as in the original. Differences from the
#' original implementation: the chi-square test is performed by
#' HardyWeinberg::HWChisq (p-values can differ marginally from the original
#' hand-rolled test), and the per-population genotype counts are computed on
#' the correct individuals (the original subset individuals with a recycled
#' per-population index).

#' @param x A genlight object containing the SNP genotypes [required].
#' @param Yates Whether to use Yates's continuity correction [default FALSE].
#' @param mono.rm Remove monomorphic loci after filtering [default FALSE].
#' @param recalc Recalculate the locus metadata statistics after filtering
#' [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @return The genlight object with excessively-heterozygous loci removed,
#' returned invisibly.

#' @author Author(s): Jesús Castrejón-Figueroa, Diana A Robledo-Ruiz (Custodian: Ching Ching Lau) -- Post
#' to \url{https://groups.google.com/d/forum/dartr}

#' @references
#' \itemize{
#' \item Robledo-Ruiz, D. A., Austin, L., Amos, J. N., Castrejon-Figueroa, J.,
#' Harley, D. K., Magrath, M. J., Sunnucks, P. & Pavlova, A. (2023).
#' Easy-to-use R functions to separate reduced-representation genomic
#' datasets into sex-linked and autosomal loci, and conduct sex assignment.
#' Molecular Ecology Resources.
#' }

#' @examples
#' \donttest{
#' if (isTRUE(getOption("dartR_fbm"))) LBP <- gl.gen2fbm(LBP)
#' filtered.gl <- gl.filter.excess.het(x = LBP, Yates = TRUE)
#' }

#' @seealso \code{\link{gl.filter.hwe}}, \code{\link{gl.report.hwe}}
#' @export

gl.filter.excess.het <- function(x,
                                 Yates = FALSE,
                                 mono.rm = FALSE,
                                 recalc = FALSE,
                                 verbose = NULL) {

  .Deprecated(
    new = "gl.filter.hwe",
    msg = paste(
      "gl.filter.excess.het() is deprecated and will be removed in a",
      "future release.\nUse gl.filter.hwe(direction = 'excess',",
      "test.type = 'ChiSquare', mult.comp.adj = TRUE,",
      "mult.comp.adj.method = 'fdr') instead."
    )
  )

  x2 <- gl.filter.hwe(
    x,
    subset = "each",
    n.pop.threshold = 1,
    test.type = "ChiSquare",
    mult.comp.adj = TRUE,
    mult.comp.adj.method = "fdr",
    alpha = 0.05,
    cc.val = ifelse(Yates, 0.5, 0),
    direction = "excess",
    min.hobs = 0.5,
    verbose = verbose
  )

  # Retain the original wrapper's optional post-processing
  hold.history <- x2@other$history
  if (mono.rm) {
    x2 <- gl.filter.monomorphs(x2, verbose = 0)
  }
  if (recalc) {
    x2 <- gl.recalc.metrics(x2, verbose = 0)
  }
  x2@other$history <- hold.history

  return(invisible(x2))
}
