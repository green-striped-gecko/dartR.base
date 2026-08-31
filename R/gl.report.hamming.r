#' @name gl.report.hamming
#' @title Reports pairwise Hamming distances between trimmed sequence tags and
#' the loci that gl.filter.hamming would remove
#' @family matched report

#' @description Hamming distance is calculated as the number of base
#' differences between two sequences. In the context of DArT trimmed
#' sequences, which are anchored to the left by the restriction enzyme
#' recognition sequence, sequences are compared starting immediately after the
#' common recognition sequence, over a fixed-length substring
#' (\code{min.length} bases), exactly as in
#' \code{\link{gl.filter.hamming}}.
#'
#' The function reports (a) the distribution of pairwise Hamming distances
#' among comparable loci, and (b) for a range of candidate thresholds, the
#' exact number of loci that \code{\link{gl.filter.hamming}} would remove.
#' The counts in (b) are obtained by running the same comparison engine, with
#' the same worst-to-best call-rate ordering, as the filter itself — they are
#' therefore exact, not an approximation from the distance distribution.

#' @param x Name of the genlight object containing the SNP data [required].
#' @param rs Number of bases to skip from the start of the TrimmedSequence
#' before extracting the comparison substring (i.e. restriction site length)
#' [default 5].
#' @param threshold Candidate maximum Hamming distance (number of mismatching
#' bases) to highlight in the plots and summary, matching the
#' \code{threshold} of \code{\link{gl.filter.hamming}}; an integer >= 0
#' [default 3].
#' @param min.length Length of the substring used for Hamming comparisons.
#' Longer sequences are truncated to it; only loci producing a substring of
#' exactly this length are compared, others are ignored (and would be retained
#' by the filter) [default 50].
#' @param max.threshold Largest threshold simulated in the loci-removed table
#' and plot [default 10].
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme User specified theme [default theme_dartR()].
#' @param plot.colors Vector with two color names for the borders and fill
#' [default c("#2171B5", "#6BAED6")].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()]
#' @param plot.file Filename (minus extension) for the RDS plot file
#' [Required for plot save]
#' @param tag.length Deprecated and ignored; distances are reported as base
#' counts over \code{min.length} bases, not proportions [default NULL].
#' @param probar Deprecated and ignored; the compiled engine makes a progress
#' bar unnecessary [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].

#' @details The function \code{\link{gl.filter.hamming}} removes one locus of
#' a pair whose substrings are within \code{threshold} mismatches, keeping
#' the locus with the better call rate. Because the filter is sequential (a
#' removed locus is not used to match further loci), the number of loci
#' removed at a given threshold cannot be derived from the pairwise distance
#' distribution alone; this function obtains it by running the filter's own
#' engine in simulation mode for each candidate threshold.
#'
#' Pairwise distances are computed in compiled code. When the number of
#' comparable locus pairs exceeds two million, a random sample of two million
#' pairs is used for the distribution (the loci-removed simulation is always
#' exact and unaffected); a note is printed when sampling occurs.
#'
#' This function requires the package \code{Rcpp} and a working C++ compiler
#' toolchain (Rtools on Windows, Xcode Command Line Tools on macOS), as the
#' comparison engine is compiled at run time (once per session).
#'
#' If plot.file is specified, plots are saved to the directory specified by
#' the user, or the global default working directory set by gl.set.wd() or to
#' the tempdir().
#'
#' Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }

#' @return Returns unaltered genlight object
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' \donttest{
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' gl.report.hamming(testset.gl[,1:100])
#' }

#' @seealso \code{\link{gl.filter.hamming}}

#' @importFrom stats sd
#' @import patchwork
#' @export

gl.report.hamming <- function(x,
                              rs = 5,
                              threshold = 3,
                              min.length = 50,
                              max.threshold = 10,
                              plot.display = TRUE,
                              plot.theme = theme_dartR(),
                              plot.colors = NULL,
                              plot.dir = NULL,
                              plot.file = NULL,
                              tag.length = NULL,
                              probar = NULL,
                              verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # DEPRECATED ARGUMENTS
    if ((!is.null(tag.length) || !is.null(probar)) && verbose >= 2) {
      cat(warn(
        "  Arguments tag.length and probar are deprecated and ignored:",
        "distances are now computed over min.length bases in compiled code\n"
      ))
    }
    if (verbose == 0) {plot.display <- FALSE}

    # SET WORKING DIRECTORY
    plot.dir <- gl.check.wd(plot.dir, verbose = 0)

    # SET COLOURS
    if (is.null(plot.colors)) {
      plot.colors <- c("#2171B5", "#6BAED6")
    } else {
      if (length(plot.colors) > 2) {
        if (verbose >= 2) {cat(warn("  More than 2 colors specified, only the first 2 are used\n"))}
        plot.colors <- plot.colors[1:2]
      }
    }

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)

    # FUNCTION SPECIFIC ERROR CHECKING

    if (!(requireNamespace("Rcpp", quietly = TRUE))) {
      stop(error(
        "Package Rcpp needed for this function to work. Please install it.\n"
      ))
    }

    if (length(x@other$loc.metrics$TrimmedSequence) == 0) {
        stop(error("Fatal Error: Data must include Trimmed Sequences\n"))
    }

    if (rs < 0 | rs > min.length) {
        stop(
            error(
                "Fatal Error: Length of restriction enzyme recognition sequence
                must be greater than zero, and less than the maximum length of a
                sequence tag; usually it is less than 9\n"
            )
        )
    }

    if (nLoc(x) == 1) {
        stop(error("Fatal Error: Data must include more than one locus\n"))
    }

    if (threshold < 0) {
      stop(error("Fatal Error: threshold must be a non-negative number of bases\n"))
    }

    if (threshold > 0 && threshold < 1) {
      stop(error(
        "Fatal Error: threshold is the maximum number of mismatching bases",
        " (e.g. 3), not a proportion, consistent with gl.filter.hamming\n"
      ))
    }

    max.threshold <- max(max.threshold, threshold)
    max.threshold <- min(max.threshold, min.length - 1)

    # DO THE JOB

    engine <- utils.hamming.engine()

    nL <- nLoc(x)
    seqs <- toupper(as.character(x@other$loc.metrics$TrimmedSequence))
    trimmed <- substr(seqs, rs + 1, min.length + rs)
    raws <- lapply(trimmed, charToRaw)
    lens <- lengths(raws)

    idx <- which(lens == min.length)
    n.short <- nL - length(idx)
    if (verbose >= 2 && n.short > 0) {
      cat(report(
        " ", n.short, "loci with a TrimmedSequence shorter than",
        min.length + rs, "bases were not compared and would be retained\n"
      ))
    }

    if (length(idx) < 2) {
      stop(error(
        "Fatal Error: fewer than two loci have a TrimmedSequence of at least",
        min.length + rs, "bases; nothing to compare\n"
      ))
    }

    # Pairwise Hamming distances among comparable loci (base counts).
    # All pairs when feasible, otherwise a random sample of max.pairs.
    if (verbose >= 2) {
        cat(report(
            "  Calculating pairwise Hamming distances between trimmed",
            "Reference sequence tags\n"
        ))
    }

    n.comp <- length(idx)
    n.pairs.total <- n.comp * (n.comp - 1) / 2
    max.pairs <- 2e6
    if (n.pairs.total <= max.pairs) {
      pairs <- t(utils::combn(idx, 2))
      sampled <- FALSE
    } else {
      i.samp <- sample(idx, max.pairs, replace = TRUE)
      j.samp <- sample(idx, max.pairs, replace = TRUE)
      ok <- i.samp != j.samp
      pairs <- cbind(i.samp[ok], j.samp[ok])
      sampled <- TRUE
      if (verbose >= 2) {
        cat(report(
          "  More than", format(max.pairs, big.mark = ","), "locus pairs;",
          "the distance distribution is based on a random sample of pairs\n"
        ))
      }
    }
    d <- engine$pairwise(raws, pairs)

    # Exact simulation of gl.filter.hamming across candidate thresholds:
    # same engine, same worst-to-best call-rate ordering as the filter.
    na.counts <- glNA(x)
    ord <- idx[order(na.counts[idx], decreasing = TRUE)]

    thresholds <- 0:max.threshold
    removed <- integer(length(thresholds))
    capped.any <- FALSE
    for (k in seq_along(thresholds)) {
      res <- engine$dedup(raws[ord], k = thresholds[k],
                          max_candidates_cap = 5000)
      removed[k] <- sum(!res$keep)
      capped.any <- capped.any || res$capped
    }
    if (capped.any && verbose >= 1) {
      cat(warn(
        "  Warning: more than 5000 candidate matches for at least one locus;",
        " removal counts at the higher thresholds may be underestimates\n"
      ))
    }

    df <- data.frame(
      Threshold = thresholds,
      Removed = removed,
      Percent.removed = round(removed * 100 / nL, 1),
      Retained = nL - removed,
      Percent.retained = round((nL - removed) * 100 / nL, 1)
    )

    # get title for plots
    if (datatype == "SNP") {
        title <- "SNP data (DArTSeq)\nPairwise Hamming distance between sequence tags"
    } else {
        title <- "Fragment P/A data (SilicoDArT)\nPairwise Hamming distance between sequence tags"
    }

    # Boxplot of pairwise distances
    p1 <-
        ggplot(data.frame(d = d), aes(y = d)) +
      geom_boxplot(color = plot.colors[1], fill = plot.colors[2]) +
      geom_hline(yintercept = threshold, color = "red", linewidth = 1) +
      coord_flip() +
      plot.theme +
      xlim(range = c(-1, 1)) +
      ylim(0, min.length) +
      ylab(" ") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      ggtitle(title)

    # Histogram of pairwise distances
    p2 <-
        ggplot(data.frame(d = d), aes(x = d)) +
      geom_histogram(binwidth = 1, color = plot.colors[1], fill = plot.colors[2]) +
      geom_vline(xintercept = threshold, color = "red", linewidth = 1) +
      coord_cartesian(xlim = c(0, min.length)) +
      xlab("Hamming distance (mismatching bases)") +
      ylab("Count") +
      plot.theme

    # Loci removed by gl.filter.hamming at each candidate threshold
    df.plot <- df
    df.plot$current <- df.plot$Threshold == threshold
    p4 <-
        ggplot(df.plot, aes(x = Threshold, y = Removed, fill = current)) +
      geom_col(color = plot.colors[1], show.legend = FALSE) +
      scale_fill_manual(values = c(`FALSE` = plot.colors[2], `TRUE` = "red")) +
      scale_x_continuous(breaks = thresholds) +
      xlab("Threshold (mismatching bases)") +
      ylab("Loci removed") +
      ggtitle("Loci that gl.filter.hamming would remove") +
      plot.theme

    if (verbose >= 3) {
      cat("    No. of loci =", nL, "\n")
      cat("    No. of individuals =", nInd(x), "\n")
      cat("    Loci compared =", n.comp, "\n")
      cat("    Minimum Hamming distance:", min(d), "bases\n")
      cat("    Maximum Hamming distance:", max(d), "bases\n")
      cat(paste0(
          "    Mean Hamming distance ",
          round(mean(d), 2),
          " +/- ",
          round(sd(d), 3),
          " SD bases\n"
      ))
      if (sampled) {
        cat("    (distance summaries from a random sample of",
            format(nrow(pairs), big.mark = ","), "pairs)\n")
      }
      cat("\n")
    }

    # PRINTING OUTPUTS
    # using package patchwork
    p3 <- (p1 / p2 / p4) + plot_layout(heights = c(1, 3, 3))
    if (plot.display) {print(p3)}
    if (verbose >= 2) {
      cat(report("  Loci removed by gl.filter.hamming at candidate thresholds",
                 "(exact simulation):\n"))
      print(df, row.names = FALSE)
    }

    if (!is.null(plot.file)) {
      tmp <- utils.plot.save(p3,
                             dir = plot.dir,
                             file = plot.file,
                             verbose = verbose)
    }

    # FLAG SCRIPT END

    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }

    # RETURN
    invisible(x)

}
