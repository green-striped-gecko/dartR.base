#' @name gl.report.bases
# PRELIMINARIES -- Set parameters --------------
#' @title Reports summary of base pair frequencies
#' @family matched reports

#' @description
#' This script calculates the frequencies of the four DNA nucleotide bases:
#' adenine (A), cytosine (C), 'guanine (G) and thymine (T), and the frequency of
#' transitions (Ts) and transversions (Tv) in a DArT genlight object.

#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param plot.display If TRUE, histograms of base composition are displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors List of two color names for the borders and fill of the
#'  plots [default c("#2171B5","#6BAED6")].
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#' @param ... Parameters passed to function \link[ggplot2]{ggsave}, 
#'  such as width and height, when the ggplot is to be saved.

#' @details 
#' The function checks first if trimmed sequences are included in the
#' locus metadata (@@other$loc.metrics$TrimmedSequence), and if so, tallies up
#' the numbers of A, T, G and C bases. Only the reference state at the SNP locus
#' is counted. Counts of transitions (Ts) and transversions (Tv) assume that
#' there is no directionality, that is C->T is the same as T->C, because the
#' reference state is arbitrary.

#' For presence/absence data (SilicoDArT), it is not possible to count
#' transversions or transitions or transversions/transitions ratio because the
#'  SNP data are not available, only a single sequence tag per locus.

#'  A color vector can be obtained with gl.select.colors() and then passed to the function
#'  with the plot.colors parameter.

#' Themes can be obtained from in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'   If a plot.file is given, the ggplot arising from this function is saved as an "RDS" 
#' binary file using saveRDS(); can be reloaded with readRDS(). A file name must be 
#' specified for the plot to be saved.

#'  If a plot directory (plot.dir) is specified, the ggplot binary is saved to that
#'  directory; otherwise to the tempdir(). 

#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' # SNP data
#'   out <- gl.report.bases(testset.gl)
#'   col <- gl.select.colors(select=c(6,1),palette=rainbow)
#'   out <- gl.report.bases(testset.gl,plot.colors=col)

#'   #' # Tag P/A data
#'   out <- gl.report.bases(testset.gs)
#' @export
#' @return The unchanged genlight object
#' 
# ----------------------
# Function
gl.report.bases <- function(x,
                            plot.display=TRUE,
                            plot.theme = theme_dartR(),
                            plot.colors = NULL,
                            plot.file=NULL,
                            plot.dir=NULL,
                            verbose = NULL,
                            ...) {
# PRELIMINARIES -- checking ----------------
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
    datatype <- utils.check.datatype(x, accept = c("genlight", "SNP", "SilicoDArT"), verbose = verbose)
    
    
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
    cat(paste("    C:", round(C, 2)), "\n")
    
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
    
    # Plot the results ----------------- 
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
      if (plot.display){print(p3)}

    # Optionally save the plot ---------------------

    if(!is.null(plot.file)){
      tmp <- utils.plot.save(p3,
                            dir=plot.dir,
                            file=plot.file,
                            verbose=verbose)
    }
    
# FLAG SCRIPT END ---------------
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
# ----------------------
    
    # RETURN
    invisible(x)
}
