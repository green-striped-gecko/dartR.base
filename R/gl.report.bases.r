#' @name gl.report.bases
# PRELIMINARIES -- Set parameters --------------
#' @title Reports summary of base pair frequencies
#'
#' @description
#' This script calculates the frequencies of the four DNA nucleotide bases:
#' adenine (A), cytosine (C), 'guanine (G) and thymine (T), and the frequency of
#' transitions (Ts) and transversions (Tv) in a DArT genlight object.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param plot.display If TRUE, histograms of base composition are displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors List of two color names for the borders and fill of the
#'  plots [default two_colors=c("#3B9AB2", "#78B7C5")].
#' @param save.type If specified, will direct the saved output to a file of this type [default NULL]
#' @param save.dir Directory in which to save the ggplot [default = working directory]
#' @param save.file Name for the ggsave file [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#'
#' @details 
#' The function checks first if trimmed sequences are included in the
#' locus metadata (@@other$loc.metrics$TrimmedSequence), and if so, tallies up
#' the numbers of A, T, G and C bases. Only the reference state at the SNP locus
#' is counted. Counts of transitions (Ts) and transversions (Tv) assume that
#' there is no directionality, that is C->T is the same as T->C, because the
#' reference state is arbitrary.
#'
#' For presence/absence data (SilicoDArT), it is not possible to count
#' transversions or transitions or transversions/transitions ratio because the
#'  SNP data are not available, only a single sequence tag per locus.
#'  
#'  A color vector can be obtained with gl.select.colors() and then passed to the function
#'  with the plot.colors parameter.
#'
#' Themes can be obtained from in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'  
#'  If the save.type parameter is set to one of 
#' "RDS", "eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", 
#' "bmp", "svg" or "wmf" (windows only), the graphics produced by the function
#' will be saved to disk. The option "RDS" saves as a binary file 
#' using saveRDS(); can be reloaded with readRDS().
#' 
#' Optional additional parameters for ggsave() can be added to the parameter list
#' (...) to govern aspects of the saved plot. Refer to ?ggsave for details.
#' 
#'  If a plot directory (save.dir) is specified, the ggplot binary is saved to that
#'  directory; otherwise to the working directory. 
#'  
#'  A file name must be specified for the plot to be saved.
#'
#' @family dartR-base
#' @export
#' @return The unchanged genlight object
#' 
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' # SNP data
#'   out <- gl.report.bases(testset.gl)
#'   out <- gl.report.bases(testset.gl,save.type="pdf",save.file="myplot")
#'   
#'   col <- gl.select.colors(select=c(6,1),palette=rainbow)
#'   out <- gl.report.bases(testset.gl,plot.colors=col)

#'   #' # Tag P/A data
#'   out <- gl.report.bases(testset.gs)
#'
# ----------------------
# Function
gl.report.bases <- function(x,
                            plot.display=TRUE,
                            plot.theme = theme_dartR(),
                            plot.colors = c("#3B9AB2", "#78B7C5"),
                            save.type=NULL,
                            save.dir=NULL,
                            save.file=NULL,
                            verbose = NULL,
                            ...) {
# PRELIMINARIES -- checking ----------------
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # # CHECK DIRECTORY FOR PLOTS
    # if(plot.save){
    #   if(is.null(save.dir)){save.dir <- tempdir()}
    # }
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    #plot.theme <- theme_dartR()
    
    if (!any(names(x@other$loc.metrics) == "TrimmedSequence")) {
        stop(error(
            "  Fatal Error: Dataset does not include variable 
            TrimmedSequence!\n"
        ))
    }
    
# DO THE JOB -- SNP data ----------------------
    
    # Count up the number of bases, and the number of each of ATGC, and other
    if (verbose >= 2) {
        cat(report("  Counting the bases\n"))
    }
    
    A <- sum(stringr::str_count(x@other$loc.metrics$TrimmedSequence, "A"))
    G <- sum(stringr::str_count(x@other$loc.metrics$TrimmedSequence, "G"))
    C <- sum(stringr::str_count(x@other$loc.metrics$TrimmedSequence, "C"))
    T <- sum(stringr::str_count(x@other$loc.metrics$TrimmedSequence, "T"))
    total <- sum(stringr::str_length(x@other$loc.metrics$TrimmedSequence))
    total.ATGC <- sum(A, G, C, T)
    if (verbose >= 2) {
        if (total != total.ATGC) {
            cat(warn("  Warning: Codes other than A, T, G and C present\n"))
        }
    }
    other <- total - total.ATGC
    other <- other * 100 / total
    A <- A * 100 / total
    G <- G * 100 / total
    T <- T * 100 / total
    C <- C * 100 / total
    
    # Calculate the fragment lengths
    mn <- mean(stringr::str_length(x@other$loc.metrics$TrimmedSequence))
    mx <- max(stringr::str_length(x@other$loc.metrics$TrimmedSequence))
    mi <- min(stringr::str_length(x@other$loc.metrics$TrimmedSequence))
    
    if (datatype == "SNP") {
        # Extract the SNPs
        matrix <- stringr::str_split_fixed(x@other$loc.metrics$SNP, ":", 2)
        state.change <- matrix[, 2]
        
        if (verbose >= 2) {
            cat(report("  Counting Transitions and Transversions\n"))
        }
        
        # Sum the transitions and transversions
        tv <-
            sum(str_count(state.change, "A>C")) + 
          sum(stringr::str_count(state.change, "C>A")) +
          sum(stringr::str_count(state.change, "G>T")) +
          sum(stringr::str_count(state.change, "T>G")) + 
          sum(stringr::str_count(state.change, "A>T")) + 
          sum(stringr::str_count(state.change, "T>A")) + 
          sum(stringr::str_count(state.change, "G>C")) + 
          sum(stringr::str_count(state.change, "C>G"))
        
        ts <-
            sum(stringr::str_count(state.change, "A>G")) + 
          sum(stringr::str_count(state.change, "G>A")) + 
          sum(stringr::str_count(state.change, "C>T")) + 
          sum(stringr::str_count(state.change, "T>C"))
        
        if (verbose >= 2) {
            if (ts + tv != length(x@other$loc.metrics$TrimmedSequence)) {
                cat(
                    warn(
                        "  Warning: Sum of transitions plus transversions does 
                        not equal number of loci.\n"
                    )
                )
            }
        }
        ts <- ts * 100 / length(x@other$loc.metrics$TrimmedSequence)
        tv <- tv * 100 / length(x@other$loc.metrics$TrimmedSequence)
        ratio <- ts / tv
    }
    
    # Printing outputs -----------
    cat(paste("  Average trimmed sequence length:",
        round(mn, digits = 1),"(",mi,"to",mx,")"),"\n")
    cat(paste(
        "  Total number of trimmed sequences:",
        length(x@other$loc.metrics$TrimmedSequence)
    ), "\n")
    cat("  Base frequencies (%)\n")
    cat(paste("    A:", round(A, 2)), "\n")
    cat(paste("    G:", round(G, 2)), "\n")
    cat(paste("    T:", round(T, 2)), "\n")
    cat(paste("    C:", round(C, 2)), "\n\n")
    
# DO THE JOB -- Tag P/A data ----------------------
    
    if (datatype == "SilicoDArT") {
        if (verbose >= 2) {
            cat(
                important(
                    "  Tag P/A data (SilicoDArT), transition/transversions 
                    cannot be calculated\n"
                )
            )
        }
        tv <- NA
        ts <- NA
    } else {
        cat(paste("  Transitions  :", round(ts, 2), "\n"))
        cat(paste("  Transversions:", round(tv, 2), "\n"))
        cat(paste("  tv/ts ratio:", round(ratio, 4), "\n\n"))
    }
    
# PLOT THE RESULTS ----------------- 
    if (plot.display | !is.null(save.type)) {
      if (datatype == "SNP") {
        title <- paste0("SNP: Base Frequencies")
      } else {
        title <- paste0("Tag P/A: Base Frequencies")
      }
      
      bases <- c("A", "C", "T", "G")
      freq <- round(c(A, C, T, G), 1)
      df <- data.frame(bases = bases, freq = freq)
      
      p1 <-
        ggplot(data = df, aes(x = bases, y = freq)) +
        geom_bar(stat="identity",color=plot.colors[1],fill=plot.colors[2]) + 
        xlab("Bases") +
        ylab("Percent Frequency") + 
        ggtitle(title) +
        plot.theme
      
      if (datatype == "SNP") {
        bases <- c("Ts", "Tv")
        freq <- round(c(ts, tv), 1)
        df2 <- data.frame(bases = bases, freq = freq)
        
        p2 <-
          ggplot(data = df2, aes(x = bases, y = freq)) +
          geom_bar(stat="identity",color=plot.colors[1],fill=plot.colors[2]) +
          xlab("Mutation Type") + 
          ylab("Percent Frequency") + 
          ggtitle(paste("SNP: Ts/Tv Rates [ratio =",round(ratio,2),"]")) +
          plot.theme
        
        p3 <- (p1 / p2)  # Using package patchwork
      } else {
        p3 <- p1
      }
      print(p3)
    }
    
    # Optionally save the plot ---------------------

    tmp <- utils.ggplotsave(p3,
                            type=save.type,
                            dir=save.dir,
                            file=save.file,
                            verbose=verbose)
# FINISH UP -------------------
    # Create return list
    if (verbose >= 2) {
        cat(
            report(
                "  Returning a list containing
[[1]] $freq -- the table of base frequencies and transition/transversion ratios;
[[2]] $plotbases -- ggplot bargraph of base frequencies;
[[3]] $plottstv -- ggplot bargraph of transitions and transversions."
            )
        )
    }
    
    out <-
        c(round(A, 2),
          round(G, 2),
          round(T, 2),
          round(C, 2),
          round(tv, 2),
          round(ts, 2))
    names(out) <- c("A", "G", "T", "C", "tv", "ts")
    
    
    # FLAG SCRIPT END 
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
# ----------------------
    
    # RETURN
    invisible(x)
}
