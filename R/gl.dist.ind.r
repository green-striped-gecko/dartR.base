#' @name gl.dist.ind
#' @title Calculates a distance matrix for individuals defined in a genlight object
#' @family distance

#' @description
#' Calculates various distances between individuals based on allele
#' frequencies or presence-absence data 

#' @param x Name of the genlight [required].
#' @param method Specify distance measure [SNP: Euclidean; P/A: Simple].
#' @param scale If TRUE, the distances are scaled to fall in the range [0,1] [default TRUE]
#' @param swap If TRUE and working with presence-absence data, then presence 
#' (no disrupting mutation) is scored as 0 and absence (presence of a disrupting 
#' mutation) is scored as 1 [default FALSE].
#' @param type Specify the type of output, dist or matrix [default dist]
#' @param plot.display If TRUE, resultant plots are displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors List of two color names for the borders and fill of the
#'  plots [default c("#2171B5","#6BAED6")].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log ; 3, progress and results summary; 5, full report
#'   [default 2 or as specified using gl.set.verbosity].
#'   
#' @details
#' The distance measure for SNP genotypes can be one of:
#' \itemize{
#'  \item Euclidean Distance [method = "Euclidean"]
#'  \item Scaled Euclidean Distance [method='Euclidean", scale=TRUE]
#'  \item Simple Mismatch Distance [method="Simple"]
#'  \item Absolute Mismatch Distance [method="Absolute"]
#'  \item Czekanowski (Manhattan) Distance [method="Manhattan"]
#'  }

#' The distance measure for Sequence Tag Presence/Absence data (binary) can be one of:
#' \itemize{
#'  \item Euclidean Distance [method = "Euclidean"]
#'  \item Scaled Euclidean Distance [method='Euclidean", scale=TRUE]
#'  \item Simple Matching Distance [method="Simple"]
#'  \item Jaccard Distance [method="Jaccard"]
#'  \item Bray-Curtis Distance [method="Bray-Curtis"]
#'  }

#' Refer to the documentation of functions in
#'   https://doi.org/10.1101/2023.03.22.533737 for algorithms
#'   and definitions.
#' 
#' @author Author(s): Custodian: Arthur Georges -- Post to #' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#'  \donttest{
#' D <- gl.dist.ind(testset.gl[1:20,], method='manhattan')
#' D <- gl.dist.ind(testset.gs[1:20,], method='Jaccard',swap=TRUE)
#' }
#' D <- gl.dist.ind(testset.gl[1:20,], method='euclidean',scale=TRUE)
#' 
#' @export
#' @return An object of class 'matrix' or dist' giving distances between individuals

gl.dist.ind <- function(x,
                        method = NULL,
                        scale = FALSE,
                        swap=FALSE,
                        type="dist",
                        plot.display = TRUE,
                        plot.theme = theme_dartR(),
                        plot.colors = NULL,
                        plot.file=NULL,
                        plot.dir=NULL,
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
    datatype <-
        utils.check.datatype(x,
                             accept = c("SNP", "SilicoDArT"),
                             verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (is.null(method) && datatype == "SNP") {
        method <- "Euclidean"
    }
    if (is.null(method) && datatype == "SilicoDArT") {
        method <- "Simple"
    }
    method <- tolower(method)
    
    if (!(
        method %in% c(
            "euclidean",
            "simple",
            "manhattan",
            "jaccard",
            "bray-curtis",
            "czekanowski",
            "absolute"
        )
    )) {
        
        if (datatype == "SNP") {
            method <- "euclidean"
            cat(warn(" Warning: Method not in the list of options, set to Euclidean Distance\n"))
        }
        if (datatype == "SilicoDArT") {
            method <- "simple"
            cat(warn(" Warning: Method not in the list of options, set to Simple Matching Distance\n"))
        }
    }
    
    # DO THE JOB
    
    if (datatype == "SNP") {
        # Calculate euclidean distance using dist 
        if (method == "euclidean") {
            if(scale==TRUE){
                dd <- utils.dist.ind.snp(x, method='euclidean',scale=TRUE,verbose=0)
                if (verbose >= 2) {
                    cat(report("  Calculating scaled Euclidean Distances between individuals\n"))
                }
            } else {
                dd <- utils.dist.ind.snp(x, method='euclidean',scale=FALSE,verbose=0)
                if (verbose >= 2) {
                    cat(report("  Calculating raw Euclidean Distances between individuals\n"))
                }
            }
        }
        
        # Calculate simple matching distance
        if (method == "simple") {
            dd <- dd <- utils.dist.ind.snp(x, method='simple',verbose=0)
            if (verbose >= 2) {
                cat(report(
                    "  Calculating simple matching distance\n"
                ))
            }
        }
        # Calculate absolute Manhattan distance
        if (method == "manhattan") {
            dd <- dd <- utils.dist.ind.snp(x, method='manhattan',verbose=0)
            if (verbose >= 2) {
                cat(report(
                    "  Calculating Manhattan distance\n"
                ))
            }
        }     
        
        # Calculate absolute Czekanowski distance
        if (method == "czekanowski") {
            dd <- dd <- utils.dist.ind.snp(x, method='czekanowski',verbose=0)
            if (verbose >= 2) {
                cat(report(
                    "  Calculating Czekanowski distance\n"
                ))
            }
        }        
        
        # Calculate absolute matching distance
        if (method == "absolute") {
            dd <- dd <- utils.dist.ind.snp(x, method='absolute',verbose=0)
            if (verbose >= 2) {
                cat(report(
                    "  Calculating absolute matching distance\n"
                ))
            }
        }        
        
        dd <- as.dist(dd)
        
        # # Revert to original order ord <- rank(pop(x)) mat <- as.matrix(dd)[ord, ord] dd <- as.dist(mat)
        mat <- as.matrix(dd)
    }
    
    if (datatype == "SilicoDArT") {
        if (method == "euclidean" && scale==FALSE) {
            if (verbose >= 2) {
                cat(report("  Calculating the Unscaled Euclidean Distances\n"))
            }
        }
        if (method == "euclidean" && scale==TRUE) {
            if (verbose >= 2) {
                cat(report("  Calculating the Scaled Euclidean Distances\n"))
            }
        }
        if (method == "simple") {
            if (verbose >= 2) {
                cat(report("  Calculating the Simple Matching Distances\n"))
            }
        }
        if (method == "jaccard") {
            if (verbose >= 2) {
                cat(report("  Calculating distances based on the Jaccard Coefficient\n"))
            }
        }
        if (method == "bray-curtis") {
            if (verbose >= 2) {
                cat(report(
                    "  Calculating the Bray-Curtis Distance\n"
                ))
            }
        }

        mat <- utils.dist.binary(x, 
                                 method = method, 
                                 swap=swap, 
                                 type="matrix", 
                                 scale = scale, 
                                 verbose = 0)
        dd <- as.dist(mat)
    }
    
    # PLOT
        if (datatype == "SNP") {
            title_plot <-
                paste0("SNP data (DArTSeq)\nInter-individual ",
                       method,
                       " distance")
        } else {
            if(method=="euclidean" && scale == TRUE){
                title_plot <-
                paste0(
                    "Presence[1]/Absence[0] data (SilicoDArT)\nInter-individual scaled ",
                    method,
                    " distance"
                )
            } else {
                if(swap==TRUE){
                    title_plot <- paste0(
                        "Presence[0]/Absence[1] data (SilicoDArT swapped)\nInter-individual ",
                        method, " distance")
                } else {
                    title_plot <- paste0(
                        "Presence[1]/Absence[0] data (SilicoDArT)\nInter-individual ",
                        method, " distance")
                }
            }
        }
        values <- NULL
        df_plot <- data.frame(values = as.vector(mat))
 
        # Boxplot
        p1 <-
            ggplot(df_plot, aes(y = values)) + 
            geom_boxplot(color = plot.colors[1], 
            fill = plot.colors[2]) + 
            coord_flip() + 
            plot.theme + 
            xlim(range = c(-1,1)) + 
            ylim(min(df_plot$values, na.rm = TRUE),
            max(df_plot$values, na.rm = TRUE)) + 
            ylab(" ") + 
            theme(axis.text.y = element_blank(), 
                  axis.ticks.y = element_blank()) + 
            ggtitle(title_plot)
        
        # Histogram
        p2 <-
            ggplot(df_plot, aes(x = values)) + 
            geom_histogram(bins = 100,
                           color = plot.colors[1],
                           fill = plot.colors[2]) + 
            xlim(min(df_plot$values, na.rm = TRUE), 
                 max(df_plot$values, na.rm = TRUE)) + 
            xlab("Distance") + 
            ylab("Count") + 
            plot.theme
        
        # PRINTING OUTPUTS
        
            # using package patchwork
            p3 <- (p1 / p2) + plot_layout(heights = c(1, 4))
            if (plot.display) {suppressWarnings(print(p3))}
            
            # Optionally save the plot ---------------------
            
            if(!is.null(plot.file)){
              tmp <- utils.plot.save(p3,
                                     dir=plot.dir,
                                     file=plot.file,
                                     verbose=verbose)
            }
    
    # SUMMARY Print out some statistics
    if (verbose >= 3) {
        cat("  Reporting inter-individual distances\n")
        cat("  Distance measure:", method, "\n")
        cat("    No. of populations =", nPop(x), "\n")
        cat("    Average no. of individuals per population =",
            round(nInd(x) / nPop(x), 1),
            "\n")
        cat("    No. of loci =", nLoc(x), "\n")
        cat("    Minimum Distance: ", round(min(dd), 2), "\n")
        cat("    Maximum Distance: ", round(max(dd), 2), "\n")
        cat("    Average Distance: ", round(mean(dd), 3), "\n\n")
    }
    
    # FLAG SCRIPT END
    
    if(type=="matrix"){
        if(verbose >= 2){
            cat(report("  Returning a square matrix\n"))
        }
        dimnames(mat) <- list(indNames(x), indNames(x))
        final <- mat
    }
    if(type!="matrix"){
        if(verbose >= 2){
            cat(report("  Returning a stat::dist object\n"))
        }
        dm <- as.matrix(dd)
        dimnames(dm) <- list(indNames(x), indNames(x))
        dd <- as.dist(dm)
        final <- dd
    }
    
    if (verbose > 0) {
            cat(report("Completed:", funname, "\n"))
    }
    return(final)
 }
