#' @name gl.dist.pop
#' @title Calculates a distance matrix for populations with SNP or Silicodart genotypes in a
#'  genlight object
#' @family distance

#' @description
#' This script calculates various distances between populations based on allele
#' frequencies (SNP genotypes) or frequency of presences in PA (SilicoDArT) data 
#'  
#' @param x Name of the genlight object [required].
#' @param as.pop Temporarily assign another locus metric as the population for
#' the purposes of deletions [default NULL].
#' @param method Specify distance measure [default euclidean].
#' @param plot.display If TRUE, resultant plots are displayed in the plot window
#' [default TRUE].
#' @param scale If TRUE and method='Euclidean', the distance will be scaled to 
#'  fall in the range [0,1] [default FALSE].
#' @param type Specify the type of output, dist or matrix [default 'dist']
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
#' For SNP data, the distance measure can be one of 'euclidean', 'fixed-diff', 'reynolds',
#' 'nei' and 'chord'. For SilicoDArT data, the distance measure can be one of 'Refer to the documentation of functions in
#'   https://doi.org/10.1101/2023.03.22.533737 for algorithms
#'   and definitions.
#'   
#' @author author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#'  \donttest{
#' # SNP genotypes
#' D <- gl.dist.pop(possums.gl, method='euclidean')
#' D <- gl.dist.pop(possums.gl, method='euclidean',scale=TRUE)
#' D <- gl.dist.pop(possums.gl, method='nei')
#' D <- gl.dist.pop(possums.gl, method='reynolds')
#' D <- gl.dist.pop(possums.gl, method='chord')
#' D <- gl.dist.pop(possums.gl, method='fixed-diff')
#' #Presence-Absence data [only 10 individuals due to speed]
#' D <- gl.dist.pop(testset.gs[1:10,], method='euclidean')
#' }
#' 
#' @export
#' @return An object of class 'dist' giving distances between populations

gl.dist.pop <- function(x,
                        as.pop=NULL,
                        method = "euclidean",
                        scale = FALSE,
                        type = "dist",
                        plot.display=TRUE,
                        plot.theme = theme_dartR(),
                        plot.colors = NULL,
                        plot.file=NULL,
                        plot.dir=NULL,
                        verbose = NULL) {

    # CHECK IF PACKAGES ARE INSTALLED
    pkg <- "reshape2"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    if(verbose==0){plot.display=FALSE}
    
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
        utils.check.datatype(x, accept = c("SNP","SilicoDArT"), verbose = verbose)
    
    if (!is(x, "dartR")) {
      class(x) <- "dartR"  
      if (verbose>2) {
        cat(warn("Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"))
        cat(warn("                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR, 
                 please use function dartR2gl\n"))
      }
    }
    
    # Population labels assigned?
    if (is.null(as.pop)) {
      if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
        if (verbose >= 2) {
          cat(
            warn(
              "  Warning: Population assignments not detected, running compliance check\n"
            )
          )
        }
        x <- gl.compliance.check(x, verbose = 0)
      }
    }
    
    # Assign the new population list if as.pop is specified 
    pop.hold <- pop(x)
    if (!is.null(as.pop)) {
      if (as.pop %in% names(x@other$ind.metrics)) {
        pop(x) <- unname(unlist(x@other$ind.metrics[as.pop]))
        if (verbose >= 2) {
          cat(report("  Temporarily assigning",as.pop,"as population\n"))
        }
      } else {
        stop(error("Fatal Error: individual metric assigned to 'pop' does not exist. Check names(gl@other$loc.metrics) and select again\n"))
      }
    }
    
    # DO THE JOB
    
    available_methods <-
        c(
            "euclidean",
            "nei",
            "reynolds",
            "chord",
            "fixed-diff",
            "simple",
            "jaccard",
            "sorensen"
        )
    
    method <- tolower(method)
    if (!(method %in% available_methods)) {
        cat(error("Fatal Error: Specified distance method is not among those 
                available (",available_methods,"), set to Euclidean.\n"))
      method <- "euclidean"

    }
    # hard.min.p <- 0.25

    nI <- nInd(x)
    nL <- nLoc(x)
    nP <- nPop(x)
    dd <- array(NA, c(nPop(x), nPop(x)))
    
    # Calculate distances 
    # if (method %in% distmethod) {
        if (verbose >= 2) {
            cat(report(paste(
                "  Calculating distances: ", method, "\n"
            )))
        }
        if(verbose >= 3){
            cat(report(
                "  Refer to https://doi.org/10.1101/2023.03.22.533737 for algorithms\n"
            ))
        }
    
    # Calculate allele frequencies for each population and locus
    f <- gl.allele.freq(x,percent=TRUE,by="popxloc",verbose=0)
    # Select only pop, locus, frequency columns
    f <- f[, c("popn", "locus", "frequency")]
    # Convert to a pop x locus matrix
    f <- reshape2::dcast(f, popn ~ locus, value.var = "frequency")
    # Reassign names to the populations, and convert from percentages to proportions
    row.names(f) <- f[, 1]
    f <- f[,-c(1)]
    p <- f / 100

# For both DArTseq and SilicoDArT
    if (method == "euclidean") {
        for (i in (1:(nP - 1))) {
            for (j in ((i + 1):nP)) {
                p_ind1 <- p[i,]
                p_ind2 <- p[j,]
                sq <- (p_ind1-p_ind2)**2
                sq <- sq[!is.na(sq)]
                L <- length(sq)
                if(scale==TRUE){
                    if(datatype=="SNP"){
                      dd[j,i] <- 0.5*sqrt(sum(sq)/L)
                    } else {
                      dd[j,i] <- sqrt(sum(sq)/L)
                    }
                } else {
                    dd[j,i] <- sqrt(sum(sq))
                }
            }
        }
    }
    # # Test code
    # x <- gl2gi(testset.gl)
    # x <- adegenet::genind2genpop(x)
    # D_check <- adegenet::dist.genpop(x,4) # Rogers D
    # hist(D_check,breaks=50)
    # D <- gl.dist.pop(testset.gl, method='euclidean',type="matrix",scale=TRUE)
    # D[upper.tri(D)] <- t(D)[upper.tri(D)]
    # hist(D/2,breaks=50)
    # #VALIDATED [with minor differences, missing handling?]

    ############################################################################
    ####### DARTSEQ
    ############################################################################
    
    if (method == "reynolds") {
        if(datatype=="SilicoDArT"){
            stop(error("Fatal Error: Reynolds Distance is not available 
                       for Silicodart presence-absence data\n"))
        }
        for (i in (1:(nP - 1))) {
            for (j in ((i + 1):nP)) {
                # Pull the loci for individuals i and j
                pind1 <- p[i,]
                pind2 <- p[j,]
                # Delete the pairwise missing
                tmp <- pind1+pind2
                pind1 <- pind1[!is.na(tmp)]
                pind2 <- pind2[!is.na(tmp)]
                # Squares
                psq <- (pind1-pind2)**2
                # Repeat for q
                qind1 <- 1-pind1
                qind2 <- 1-pind2
                qsq <- (qind1-qind2)**2
                # Cross products
                p12 <- pind1*pind2
                q12 <- qind1*qind2
                # Non-missing loci
                #L <- length(psq)
                
                #dd[j,i] <- sqrt(sum(psq+qsq)/(2*sum(1-p12-q12)))
                dd[j,i] <- -log(1-sqrt(sum(psq+qsq)/(2*sum(1-p12-q12))))
                #dd[j,1] <- sqrt(sum(psq)/(sum(1-p12-q12)))
            }
        }
    }
    # # Test code
    # x <- gl2gi(testset.gl)
    # x <- adegenet::genind2genpop(x)
    # D_check <- adegenet::dist.genpop(x,3) # Reynolds in common use
    # D_check <- -log(1-D_check) # Proportional to divergence time
    # hist(D_check,breaks=50)
    # D <- gl.dist.pop(testset.gl, method='reynolds',type='matrix',scale=TRUE)
    # D[upper.tri(D)] <- t(D)[upper.tri(D)]
    # hist(D,breaks=50)
    # #VALIDATED [with minor difference, missing handling?]
    
    if (method == "nei") {
        if(datatype=="SilicoDArT"){
            stop(error("Fatal Error: Nei Standard Distance is not available
                       for Silicodart presence-absence data\n"))
        }
        for (i in (1:(nP - 1))) {
            for (j in ((i + 1):nP)) {
                # Pull the loci for individuals i and j
                prow1 <- p[i,]
                prow2 <- p[j,]
                # Delete the pairwise missing
                tmp <- prow1+prow2
                prow1 <- prow1[!is.na(tmp)]
                prow2 <- prow2[!is.na(tmp)]
                # Squares
                p1sq <- prow1*prow1
                p2sq <- prow2*prow2
                # Repeat for q=1-p
                qrow1 <- 1-prow1
                qrow2 <- 1-prow2
                q1sq <- qrow1*qrow1
                q2sq <- qrow2*qrow2
                # Cross products
                p12 <- prow1*prow2
                q12 <- qrow1*qrow2
                # Number of non-missing loci
                L <- length(p12)

                dd[j,i] <- -log(sum(p12+q12)/(sqrt(sum(p1sq+q1sq))*sqrt(sum(p2sq+q2sq))))
            }
        }
    }
    # # Test code
    # x <- gl2gi(testset.gl)
    # x <- adegenet::genind2genpop(x)
    # D_check <- adegenet::dist.genpop(x,1) 
    # hist(D_check,breaks=50)
    # D <- gl.dist.pop(testset.gl, method='nei',type='matrix',scale=TRUE)
    # hist(D,breaks=50)
    # #VALIDATED [with minor difference, missing handling?]
    
    if (method == "chord") {
        if(datatype=="SilicoDArT"){
            stop(error("Fatal Error: Czfordi-Edwards Chord Distance is not available
                       for Silicodart presence-absence data\n"))
        }
        for (i in (1:(nP - 1))) {
            for (j in ((i + 1):nP)) {
                # Pull the loci for individuals i and j
                prow1 <- p[i,]
                prow2 <- p[j,]
                # Delete the pairwise missing
                tmp <- prow1+prow2
                prow1 <- prow1[!is.na(tmp)]
                prow2 <- prow2[!is.na(tmp)]
                # create proportions for allele 2
                qrow1 <- 1-prow1
                qrow2 <- 1-prow2
                # Cross products
                p12 <- prow1*prow2
                q12 <- qrow1*qrow2
                # Non-missing Loci
                L <- length(p12)

                dd[j,i] <- (2/pi)*sqrt(2*(1 - (sum(sqrt(p12))/L + sum(sqrt(q12)/L))))
            }
        }
    }
    # # Test code
    # x <- gl2gi(testset.gl)
    # x <- adegenet::genind2genpop(x)
    # D_check <- adegenet::dist.genpop(x,2) # Angular or Edwards?
    # #D_check <- -log(1-D_check) # Proportional to divergence time
    # hist(D_check,breaks=50)
    # D <- gl.dist.pop(testset.gl, method='chord',type='matrix',scale=TRUE)
    # D[upper.tri(D)] <- t(D)[upper.tri(D)]
    # hist(D,breaks=50)
    # #VALIDATED [with minor difference, missing handling?]

    if (method == "fixed-diff") {
        dd <- gl.fixed.diff(x, verbose = 0)[[3]]/100
        if (verbose >= 2) {
            cat(report("  Calculating proportion of fixed differences\n"))
            cat(
                warn(
                    "Note: this distance may be non-metric, and so should be considered a dissimilarity measure\n"
                )
            )
        }
    }

    # # Revert to original order ord <- rank(popNames(x)) mat <- as.matrix(dd)[ord, ord] dd <- as.dist(mat)
    
    if(method != "fixed-diff") {
    dimnames(dd) <- list(popNames(x), popNames(x))
    }
    
    ############################################################################
    ####### SILICODART
    ############################################################################
    
    if (method == "simple") {
      if(datatype=="SNP"){
        stop(error("Fatal Error: Simple Matching Distance is not available
                       for SNP data\n"))
      }
      tmp <- gl.dist.ind(x,method="simple",verbose=0)
      dd <- utils.collapse.matrix(D=tmp,x=x,verbose=0)
    }
    
    if (method == "jaccard") {
      if(datatype=="SNP"){
        stop(error("Fatal Error: Jaccard Distance is not available
                       for SNP data\n"))
      }
      tmp <- gl.dist.ind(x,method="jaccard",verbose=0)
      dd <- utils.collapse.matrix(D=tmp,x=x,verbose=0)
    }
    
    if (method == "sorensen") {
      if(datatype=="SNP"){
        stop(error("Fatal Error: Sorensen (=Dice) Distance is not available
                       for SNP data\n"))
      }
      tmp <- gl.dist.ind(x,method="sorensen",verbose=0)
      dd <- utils.collapse.matrix(D=tmp,x=x,verbose=0)
    }
    
    
    ############################################################################
    ####### PLOT RESULTS
    ############################################################################

    # PLOT Plot Box-Whisker plot
    
    if (plot.display) {
        if (datatype == "SNP") {
            title_plot <- paste0("SNP data\nUsing ", method, " distance")
        } else {
            title_plot <-
                paste0("Tag P/A data (SilicoDArT)\nUsing ",
                       method,
                       " distance")
        }
        values <- NULL
        val <- as.vector(dd)
        val <- val[!is.na(val)]
        df_plot <- data.frame(values = val)
        
        # Boxplot
        p1 <- ggplot(df_plot, aes(y = values)) +
        geom_boxplot(color = plot.colors[1], fill = plot.colors[2]) +
        coord_flip()  +
        plot.theme  +
        xlim(range = c(-1,1)) + 
        ylim(min(df_plot$values, na.rm = TRUE),max(df_plot$values, na.rm = TRUE)) + 
        ylab(" ") + 
        theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) + 
        ggtitle(title_plot)
        
        # Histogram
        p2 <- ggplot(df_plot, aes(x = values)) +
        geom_histogram(bins = 20,color = plot.colors[1], fill = plot.colors[2]) +
        xlim(min(df_plot$values, na.rm = TRUE), max(df_plot$values, na.rm = TRUE)) +
        xlab("Distance Metric") +
        ylab("Count") +
        plot.theme
    }
    
    # SUMMARY Print out some statistics
    if (verbose >= 3) {
        cat("  Reporting inter-population distances\n")
        cat("  Distance measure:", method, "\n")
        cat("    No. of populations =", nPop(x), "\n")
        cat("    Average no. of individuals per population =",
            round(nInd(x) / nPop(x),1),
            "\n")
        cat("    No. of loci =", nLoc(x), "\n")
        cat("    Minimum Distance: ", round(min(dd,na.rm=TRUE), 2), "\n")
        cat("    Maximum Distance: ", round(max(dd,na.rm=TRUE), 2), "\n")
        cat("    Average Distance: ", round(mean(dd,na.rm=TRUE), 3), "\n")
    }
    
    # PRINTING OUTPUTS
    
        # using package patchwork
        
        if (plot.display) {
          p3 <- (p1 / p2) + plot_layout(heights = c(1, 4))
          suppressWarnings(print(p3))}
        
        # Optionally save the plot ---------------------
        
        if(!is.null(plot.file)){
          tmp <- utils.plot.save(p3,
                                 dir=plot.dir,
                                 file=plot.file,
                                 verbose=verbose)
        }
    
    # Reassign the initial population list if as.pop is specified
    
    if (!is.null(as.pop)) {
      pop(x) <- pop.hold
      if (verbose >= 2) {
        cat(report("  Restoring population assignments to initial state\n"))
      }
    }
        
    if(type=="dist"){
    dd <- as.dist(dd)
        if(verbose >= 2){cat(report("  Returning a stats::dist object\n"))}
    } else {
        dd <- as.matrix(dd)
        if(verbose >= 2){cat(report("  Returning a square matrix object\n"))}
    }
   
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
   
    return(dd)
}
