#' @name gl.filter.hamming
#' @title Filters loci based on pairwise Hamming distance between sequence tags
#' @family matched filter

#' @description
#' Hamming distance is calculated as the number of base differences between two
#' sequences which can be expressed as a count or a proportion. Typically, it is
#' calculated between two sequences of equal length. In the context of DArT
#' trimmed sequences, which differ in length but which are anchored to the left
#' by the restriction enzyme recognition sequence, it is sensible to compare the
#' two trimmed sequences starting from immediately after the common recognition
#' sequence and terminating at the last base of the shorter sequence.
#' 
#' @param x Name of the genlight object containing the SNP data [required].
#' @param threshold A threshold Hamming distance for filtering loci
#' [default threshold 0.2].
#' @param rs Number of bases in the restriction enzyme recognition sequence
#' [default 5].
#' @param tag.length Typical length of the sequence tags [default 69].
#' @param plot.display If TRUE, histograms are displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors List of two color names for the borders and fill of the
#'  plots [default c("#2171B5", "#6BAED6")].
#' @param plot.dir Directory in which to save files [default = working directory]
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param pb If TRUE, a progress bar will be displayed [default FALSE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log ; 3, progress and results summary; 5, full report
#'   [default 2, unless specified using gl.set.verbosity].
#'   
#' @details
#' Hamming distance can be computed
#' by exploiting the fact that the dot product of two binary vectors x and (1-y)
#' counts the corresponding elements that are different between x and y.
#' This approach can also be used for vectors that contain more than two 
#' possible values at each position (e.g. A, C, T or G).

#' If a pair of DNA sequences are of differing length, the longer is truncated.

#' The algorithm is that of Johann de Jong
#'\url{https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/}
#' as implemented in \code{\link{utils.hamming}}.

#' Only one of two loci are retained if their Hamming distance is less that a 
#' specified
#' percentage. 5 base differences out of 100 bases is a 20% Hamming distance.


#' @author Custodian: Arthur Georges -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' # SNP data
#' test <- platypus.gl
#' test <- gl.subsample.loci(platypus.gl,n=50)
#' result <- gl.filter.hamming(test, threshold=0.6, verbose=3)

#' @import patchwork
#' @export
#' @return A genlight object filtered on Hamming distance.

gl.filter.hamming <- function(x,
                              threshold = 0.2,
                              rs = 5,
                              tag.length = 69,
                              plot.display=TRUE,
                              plot.theme = theme_dartR(),
                              plot.colors = NULL,
                              plot.file=NULL,
                              plot.dir=NULL,
                              pb = FALSE,
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
    if (threshold < 0 || threshold > 1) {
        cat(
            warn(
                "  Warning: Parameter 'threshold' must be an integer between 0 
                and 1, set to 0.2\n"
            )
        )
        threshold <- 0.2
    }
    
    if (length(x@other$loc.metrics$TrimmedSequence) == 0) {
        stop(error("Fatal Error: Data must include Trimmed Sequences\n"))
    }
    
    if (nLoc(x) == 1) {
        stop(error("Fatal Error: Data must include more than one locus\n"))
    }
    
    
    # DO THE JOB
    
    n0 <- nLoc(x)
    
    if (verbose >= 3) {
        cat(
            report(
                "  Note: Hamming distance ranges from zero (sequence identity)
                to 1 (no bases shared at any position)\n"
            )
        )
        cat(
            report(
                "  Note: Calculating pairwise Hamming distances between trimmed 
                reference sequence tags\n"
            )
        )
    }
    
    x@other$loc.metrics$TrimmedSequence <-
        as.character(x@other$loc.metrics$TrimmedSequence)
    
    count <- 0
    nL <- nLoc(x)
    index <- rep(TRUE, (nL - 1))
    d <- rep(NA, (((nL - 1) * nL) / 2))
    if (pb) {
        pbar <-
            txtProgressBar(
                min = 0,
                max = 1,
                style = 3,
                initial = 0,
                label = "Working ....\n"
            )
        getTxtProgressBar(pbar)
    }
    if (verbose >= 2) {
        cat(report(
            "  Filtering loci with a Hamming Distance of less than",
            threshold,
            "\n"
        ))
    }
    for (i in 1:(nL - 1)) {
        s1 <- x@other$loc.metrics$TrimmedSequence[i]
        for (j in ((i + 1):nL)) {
            count <- count + 1
            s2 <- x@other$loc.metrics$TrimmedSequence[j]
            d[count] <- utils.hamming(s1, s2, r = rs)
            if (d[count] <= threshold) {
                index[i] <- FALSE
                if (verbose >= 3) {
                    cat(
                        " Deleting:",
                        locNames(x)[i],
                        locNames(x)[j],
                        "\n"
                    )
                }
                break
            }
        }
        if (pb) {
            setTxtProgressBar(pbar, i / (nL - 1))
        }
    }
    d <- d[!is.na(d)]
    
      x2 <- x[, (index)]
      x2@other$loc.metrics <- x@other$loc.metrics[(index), ]
    
    # PLOT HISTOGRAMS, BEFORE AFTER
    if (plot.display) {
        plotvar <- d
        # min <- min(plotvar,threshold,na.rm=TRUE) min <- trunc(min*100)/100
        max <- max(plotvar, threshold, na.rm = TRUE)
        max <- ceiling(max * 10) / 10
        if (datatype == "SNP") {
            xlabel <- "Pre-filter SNP Hamming Distance"
        } else {
            xlabel <- "Pre-filter P/A Hamming Distance"
        }
        p1 <-
            ggplot(data.frame(plotvar), aes(x = plotvar)) + 
            geom_histogram(bins = 100,
                           color = plot.colors[1],
                           fill = plot.colors[2]) + 
            coord_cartesian(xlim = c(0, max)) +
            geom_vline(xintercept = threshold,
                       color = "red",
                       size = 1) + 
            xlab(xlabel) + 
            ylab("Count") + 
            plot.theme
        
        # if (datatype=='SilicoDArT'){ rdepth <-
        #x2@other$loc.metrics$AvgReadDepth } else if 
        #(datatype=='SNP'){ rdepth <-
        # x2@other$loc.metrics$rdepth }
        plotvar <- d[d >= threshold]
        # min <- min(plotvar,threshold) min <- trunc(min*100)/100 max <- 
        #max(plotvar,threshold,na.rm=TRUE) max <- ceiling(max*10)/10
        if (datatype == "SNP") {
            xlabel <- "Post-filter SNP Hamming Distance"
        } else {
            xlabel <- "Post-filter P/A Hamming Distance"
        }
        p2 <-
            ggplot(data.frame(plotvar), aes(x = plotvar)) +
            geom_histogram(bins = 100,
                           color = plot.colors[1],
                           fill = plot.colors[2]) + 
            coord_cartesian(xlim = c(0, max)) + 
            geom_vline(xintercept = threshold,color = "red", size = 1) + 
            xlab(xlabel) +
            ylab("Count") + 
            plot.theme
        
        p3 <- (p1 / p2) + plot_layout(heights = c(1, 1))
        print(p3)
    }
      
      # Optionally save the plot ---------------------
      
      if(!is.null(plot.file)){
        tmp <- utils.plot.save(p3,
                               dir=plot.dir,
                               file=plot.file,
                               verbose=verbose)
      }
    
    # REPORT A SUMMARY
    if (verbose >= 3) {
        cat("\n  Summary of filtered dataset\n")
        cat(paste("    Initial No. of loci:", n0, "\n"))
        cat(paste(
            "    Hamming d >",
            threshold,
            "=",
            round(threshold * tag.length, 0),
            "bp\n"
        ))
        cat(paste("    Loci deleted", (n0 - nLoc(x2)), "\n"))
        cat(paste("    Final No. of loci:", nLoc(x2), "\n"))
        cat(paste("    No. of individuals:", nInd(x2), "\n"))
        cat(paste("    No. of populations: ", length(levels(
            factor(pop(x2))
        )), "\n"))
    }
    
    # ADD TO HISTORY
    nh <- length(x2@other$history)
    x2@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    return(x2)
}
