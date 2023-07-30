#' @name gl.report.hamming
#' @title Calculates the pairwise Hamming distance between DArT trimmed DNA
#' sequences
#' @family matched report

#' @description Hamming distance is calculated as the number of base differences
#' between two sequences which can be expressed as a count or a proportion.
#' Typically, it is calculated between two sequences of equal length. In the
#' context of DArT trimmed sequences, which differ in length but which are
#' anchored to the left by the restriction enzyme recognition sequence, it is
#' sensible to compare the two trimmed sequences starting from immediately after
#' the common recognition sequence and terminating at the last base of the
#' shorter sequence.

#' @param x Name of the genlight object containing the SNP data [required].
#' @param rs Number of bases in the restriction enzyme recognition sequence
#' [default 5].
#' @param threshold Minimum acceptable base pair difference for display on the
#' boxplot and histogram [default 3].
#' @param tag.length Typical length of the sequence tags [default 69].
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme User specified theme [default theme_dartR()].
#' @param plot.colors Vector with two color names for the borders and fill
#' [default c("#2171B5", "#6BAED6")].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]
#' @param plot.file Filename (minus extension) for the RDS plot file [Required for plot save]
#' @param probar If TRUE, a progress bar is displayed during run [defalut FALSE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].

#' @details The function \code{\link{gl.filter.hamming}} will filter out one of
#' two loci if their Hamming distance is less than a specified percentage

#' Hamming distance can be computed by exploiting the fact that the dot product
#' of two binary vectors x and (1-y) counts the corresponding elements that are
#' different between x and y. This approach can also be used for vectors that
#' contain more than two possible values at each position (e.g. A, C, T or G).

#' If a pair of DNA sequences are of differing length, the longer is truncated.

#' The algorithm is that of Johann de Jong
#' \url{https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/}
#' as implemented in \code{\link{utils.hamming}}

#'   If plot.file is specified, plots are saved to the directory specified by the user, or the global
#'   default working directory set by gl.set.wd() or to the tempdir().

#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }

#' @return Returns unaltered genlight object
#' @author Custodian: Arthur Georges -- Post to 
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#'  \donttest{
#' gl.report.hamming(testset.gl[,1:100])
#' gl.report.hamming(testset.gs[,1:100])
#' }

#' #' # SNP data
#' test <- platypus.gl
#' test <- gl.subsample.loci(platypus.gl,n=50)
#' result <- gl.report.hamming(test, verbose=3)
#' result <- gl.report.hamming(test, plot.file="ttest", verbose=3)

#' @seealso \code{\link{gl.filter.hamming}}

#' @importFrom stats sd
#' @import patchwork
#' @export

gl.report.hamming <- function(x,
                              rs = 5,
                              threshold = 3,
                              tag.length = 69,
                              plot.display=TRUE,
                              plot.theme = theme_dartR(),
                              plot.colors = NULL,
                              plot.dir=NULL,
                              plot.file=NULL,
                              probar = FALSE,
                              verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    if(verbose==0){plot.display <- FALSE}
    
    # SET WORKING DIRECTORY
    plot.dir <- gl.check.wd(plot.dir,verbose=0)
	
    # SET COLOURS
    if(is.null(plot.colors)){
      plot.colors <- c("#2171B5", "#6BAED6")
    } else {
      if(length(plot.colors) > 2){
        if(verbose >= 2){cat(warn("  More than 2 colors specified, only the first 2 are used\n"))}
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
    
    if (length(x@other$loc.metrics$TrimmedSequence) == 0) {
        stop(error("Fatal Error: Data must include Trimmed Sequences\n"))
    }
    
    if (rs < 0 | rs > tag.length) {
        stop(
            error(
                "Fatal Error: Length of restriction enzyme recognition sequence
                must be greater than zero, and less that the maximum length of a
                sequence tag; usually it is less than 9\n"
            )
        )
    }
    
    if (nLoc(x) == 1) {
        stop(error("Fatal Error: Data must include more than one locus\n"))
    }
    
    # DO THE JOB
    
    s <- as.character(x@other$loc.metrics$TrimmedSequence)
    tld <- threshold / (tag.length - rs)
    
    if (probar) {
        pb <-
            txtProgressBar(
                min = 0,
                max = 1,
                style = 3,
                initial = 0,
                label = "Working ...."
            )
        getTxtProgressBar(pb)
    }
    
    if (verbose >= 3) {
        cat(
            report(
                "  Hamming distance ranges from zero (sequence identity) to 1 
                (no bases shared at any position)\n"
            )
        )
    }
    if (verbose >= 2) {
        cat(
            report(
                "  Calculating pairwise Hamming distances between trimmed 
                Reference sequence tags\n"
            )
        )
    }
    
    count <-0
    nL <- nLoc(x)
    d <- rep(NA, (((nL - 1) * nL) / 2))
    
    for (i in 1:(nL - 1)) {
        for (j in ((i + 1):nL)) {
            count <- count + 1
            d[count] <- utils.hamming(s[i], s[j], r = rs)
        }
        if (probar) {
            setTxtProgressBar(pb, i / (nL - 1))
        }
    }
    
    # get title for plots
    if (datatype == "SNP") {
        title <-
            paste0("SNP data (DArTSeq)\nPairwise Hamming Distance between 
                   sequence tags")
    } else {
        title <-
            paste0(
                "Fragment P/A data (SilicoDArT)\nPairwise Hamming Distance 
                between sequence tags"
            )
    }
    
    if (verbose >= 2) {
        if (plot.display) {
            cat(
                report(
                    "  Plotting boxplot and histogram of Hamming distance, 
                    showing a threshold of",
                    threshold,
                    "bp [HD",
                    round(tld, 2),
                    "]\n"
                )
            )
        }
    }
    
    # Boxplot
    p1 <-
        ggplot(as.data.frame(d), aes(y = d)) +
      geom_boxplot(color = plot.colors[1], fill = plot.colors[2]) + 
      geom_hline(yintercept = tld,color = "red", size = 1) + 
      coord_flip() + 
      plot.theme + 
      xlim(range = c(-1, 1)) + 
      ylim(0, 1) +
      ylab(" ") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
      ggtitle(title)
    
    # Histogram
    p2 <-
        ggplot(as.data.frame(d), aes(x = d)) + 
      geom_histogram(bins = 50, color = plot.colors[1],fill = plot.colors[2]) +
      geom_vline(xintercept = tld,color = "red",size = 1) + 
      coord_cartesian(xlim = c(0, 1)) +
      xlab("Hamming distance") +
      ylab("Count") + 
      annotate(geom = "text", 
               x = tld + 0.2, 
               y = max(graphics::hist(d, breaks = seq(0, 1, by = 1 / 50),
                                      plot = FALSE)$counts) * 0.75,
               label = paste("Threshold of\n", threshold,
                             "bp [HD", round(tld, 2), "]")) + 
      plot.theme
    
    cat("    No. of loci =", nLoc(x), "\n")
    cat("    No. of individuals =", nInd(x), "\n")
    cat("    Minimum Hamming distance: ", round(min(d), 2), "\n")
    cat("    Maximum Hamming distance: ", round(max(d), 2), "\n")
    cat(paste0(
        "    Mean Hamming Distance ",
        round(mean(d), 2),
        "+/-",
        round(sd(d), 3),
        " SD\n"
    ))
    n.outliers <- sum(d <= (threshold / (tag.length - rs)))
    cat(
        "    No. of pairs with Hamming Distance less than or equal to",
        threshold,
        "base pairs:",
        n.outliers,
        "\n\n"
    )
    
    # Determine the loss of loci for a given threshold using quantiles
    nl <- nLoc(x)
    quantile_res <- quantile(d, probs = seq(0, 1, 1 / 20),type=1)
    retained <- unlist(lapply(quantile_res, function(y) {
        res <- sum(d >= y)
    }))
    pc.retained <- round(retained * 100 / ((((
        nL - 1
    ) * nL) / 2)), 1)
    filtered <- ((((nL - 1) * nL) / 2)) - retained
    pc.filtered <- 100 - pc.retained
    df <-
        data.frame(as.numeric(sub("%", "", names(quantile_res))),
                   quantile_res,
                   retained,
                   pc.retained,
                   filtered,
                   pc.filtered)
    colnames(df) <-
        c("Quantile",
          "Threshold",
          "Retained",
          "Percent",
          "Filtered",
          "Percent")
    df <- df[order(-df$Quantile), ]
    df$Quantile <- paste0(df$Quantile, "%")
    rownames(df) <- NULL
    
    # PRINTING OUTPUTS
    # using package patchwork
    p3 <- (p1 / p2) + plot_layout(heights = c(1, 4))
    if (plot.display) {print(p3)}
    print(df)
    
    if(!is.null(plot.file)){
      tmp <- utils.plot.save(p3,
                             dir=plot.dir,
                             file=plot.file,
                             verbose=verbose)
    }
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    invisible(x)
    
}
