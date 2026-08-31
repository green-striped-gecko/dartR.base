#' @name gl.report.excess.het
#' @title Report loci with excess of heterozygosity [deprecated]
#' @family unmatched report

#' @description
#' \strong{Deprecated.} Use
#' \code{gl.report.hwe(direction = 'excess', method_sig = 'ChiSquare',
#' multi_comp = TRUE, multi_comp_method = 'fdr')} instead — see
#' \code{\link{gl.report.hwe}}.
#'
#' This wrapper reproduces the workflow of Robledo-Ruiz et al. (2023) for
#' identifying loci with a significant excess of heterozygotes (candidate
#' sex-linked or paralogous loci) by calling gl.report.hwe with the
#' equivalent settings, and will be removed in a future release.

#' @details
#' The published workflow is reproduced via gl.report.hwe(direction =
#' 'excess', min.hobs = 0.5, ...): only loci with observed heterozygosity of
#' at least 0.5 are tested, and the false-discovery-rate adjustment is
#' confined to that screened set, as in the original. Differences from the
#' original implementation: the plots are gl.report.hwe's ternary plots
#' rather than the original before/after scatter plots, and the chi-square
#' test is performed by HardyWeinberg::HWChisq (p-values can differ
#' marginally from the original hand-rolled test). The parameters
#' plot.theme, plot.colors, plot.file and plot.dir are accepted for backward
#' compatibility but ignored.

#' @param x Name of the genlight object containing the SNP data [required].
#' @param Yates Boolean for Yates's continuity correction. [default FALSE]
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme Ignored (retained for backward compatibility).
#' @param plot.colors Ignored (retained for backward compatibility).
#' @param plot.dir Ignored (retained for backward compatibility).
#' @param plot.file Ignored (retained for backward compatibility).
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @return A list with two elements: 'results.table' (the gl.report.hwe
#' table for loci with significant heterozygote excess) and 'removed.loci'
#' (a vector of the names of those loci), returned invisibly.

#' @author Jesús Castrejón-Figueroa, Diana A Robledo-Ruiz (Custodian: Ching Ching Lau) -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

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
#' filtered.table <- gl.report.excess.het(x = LBP, Yates = TRUE)
#' }

#' @seealso \code{\link{gl.report.hwe}}, \code{\link{gl.filter.hwe}}
#' @export

gl.report.excess.het <- function(x,
                                 Yates = FALSE,
                                 plot.display = TRUE,
                                 plot.theme = NULL,
                                 plot.colors = NULL,
                                 plot.file = NULL,
                                 plot.dir = NULL,
                                 verbose = NULL) {

  .Deprecated(
    new = "gl.report.hwe",
    msg = paste(
      "gl.report.excess.het() is deprecated and will be removed in a",
      "future release.\nUse gl.report.hwe(direction = 'excess',",
      "method_sig = 'ChiSquare', multi_comp = TRUE,",
      "multi_comp_method = 'fdr') instead."
    )
  )

  df <- gl.report.hwe(
    x,
    subset = "each",
    method_sig = "ChiSquare",
    multi_comp = TRUE,
    multi_comp_method = "fdr",
    alpha_val = 0.05,
    cc_val = ifelse(Yates, 0.5, 0),
    direction = "excess",
    min.hobs = 0.5,
    sig_only = TRUE,
    plot.out = plot.display,
    verbose = verbose
  )

  return(invisible(list(
    results.table = as.data.frame(df),
    removed.loci = unique(as.character(df$Locus))
  )))
}
