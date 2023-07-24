#' @name gl.smearplot
#' @title Smear plot of SNP or presence/absence (SilicoDArT) data
#' @family graphics

#' @description
#' Each locus is color coded for scores of 0, 1, 2 and NA for SNP data and 0, 1
#' and NA for presence/absence (SilicoDArT) data. Individual labels can be added
#' and individuals can be grouped by population.

#' Plot may become cluttered if ind.labels If there are too many individuals, 
#' it is best to use ind.labels.size = 0.

#' @param x Name of the genlight object containing the SNP or presence/absence
#'  (SilicoDArT) data [required].
#' @param ind.labels If TRUE, individuals are labelled with indNames(x) [default FALSE].
#' @param group.pop If ind.labels is TRUE, group by population [default TRUE].
#' @param ind.labels.size Size of the individual labels [default 10].
#' @param het.only If TRUE, show only the heterozygous state [default FALSE]
#' @param plot.display If TRUE, histograms of base composition are displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors List of two color names for the borders and fill of the
#'  plots [default c("#2171B5", "#6BAED6")].
#' @param plot.dir Directory in which to save files [default = working directory]
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param legend Position of the legend: “left”, “top”, “right”, “bottom” or
#'  'none' [default = 'bottom'].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity]
#' 
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' gl.smearplot(testset.gl,ind.labels=FALSE)
#' gl.smearplot(testset.gs[1:10,],ind.labels=TRUE)

#' @export
#' @return Returns unaltered genlight object

gl.smearplot <- function(x,
                        plot.display=TRUE,
                        ind.labels = FALSE,
                        ind.labels.size = 10,
                        group.pop = FALSE, 
                        plot.theme = theme_dartR(),
                        plot.colors = NULL,
                        plot.file=NULL,
                        plot.dir=NULL,
                        het.only=FALSE,
                        legend = "bottom",
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
    
    # SET WORKING DIRECTORY
    plot.dir <- gl.check.wd(plot.dir,verbose=0)
    
    # SET COLOURS
    if(is.null(plot.colors)){
      plot.colors <- gl.select.colors(library="brewer",palette="Blues",select=c(7,5))
    }
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose) 
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)
    
    # SET COLORS
    
    if(is.null(plot.colors)){
      plot.colors <- gl.select.colors(library="brewer",palette="spectral",select=c(1,8,11,6),ncolors=11)
        #plot.colors <- gl.select.colors(library="baseR",palette="topo.colors",select=c(1,5,3,9))
        #plot.colors <- gl.select.colors(library="baseR",palette="rainbow",select=c(1,4,7,3))
        #plot.colors <- c("#a6cee3","#1f78b4","#b2df8a","#dddddd") # = default for plot()
        #plot.colors <- c("#1b9e77","#d95f02","#7570b3","#dddddd") # = default for plot()
    } 
    
    if (het.only) {
      #plot.colors <- c("#dddddd", "#ff0000", "#dddddd","#dddddd" )
      plot.colors <- gl.select.colors(library="brewer",select=c(5,9,5,6))
    }
    
    # SET IND LABELS
    
    n10 <- nchar(as.character(nInd(x)))
    lzs <- paste0("%0",as.character(n10),"d")
    if(ind.labels == TRUE){
      individuals <- paste0(sprintf(lzs,1:nInd(x)),"_",indNames(x))
    } else {
      #individuals <- paste0(sprintf(lzs,1:nInd(x)))
      individuals <- "-"
    }
    
    # DO THE JOB
    
    X_temp <- as.data.frame(as.matrix(x))
    colnames(X_temp) <- 1:nLoc(x)
    X_temp$id <- individuals
    X_temp$pop <- pop(x)
    
    X <- reshape2::melt(X_temp, id.vars = c("pop", "id"))
    X$value <- as.character(X$value)
    X$value <- ifelse(X$value=="NA", NA, X$value)
    colnames(X) <- c("pop", "id", "locus", "genotype")
    loc_labels <- pretty(1:nLoc(x), 5)
    id_labels <- pretty(1:nInd(x), 5)
    
    locus <- id <- genotype <- NA
    
    if (datatype == "SilicoDArT") {
        p3 <-
            ggplot(X, aes(
                x = locus,
                y = id,
                fill = genotype
            )) + geom_raster() + scale_fill_discrete(
                type = plot.colors[c(1, 3)],
                na.value = plot.colors[4],
                name = "Genotype",
                labels = c("0", "1")
            ) + theme_dartR() + theme(
                legend.position = legend,
                axis.text.y = element_text(size = ind.labels.size)
            ) +
            scale_x_discrete(
                breaks = loc_labels,
                labels = as.character(loc_labels),
                name = "Loci"
            ) + 
            ylab("Individuals")
    }
    
    if (datatype == "SNP") {
        p3 <-
            ggplot(X, aes(
                x = locus,
                y = id,
                fill = genotype
            )) + geom_raster() + 
                scale_fill_discrete(
                type = plot.colors,
                na.value = plot.colors[4],
                name = "Genotype",
                labels = c("0", "1", "2")
            ) + theme_dartR() + theme(
                legend.position = legend,
                axis.text.y = element_text(size = ind.labels.size)
            ) +
            scale_x_discrete(
                breaks = loc_labels,
                labels = as.character(loc_labels),
                name = "Loci",
                position="bottom"
            ) + 
        ylab("Individuals")
    }
    
    if (ind.labels==TRUE & group.pop == TRUE) {
        p3 <- p3 + facet_wrap(~ pop,
                              ncol = 1,
                              dir = "v",
                              scales = "free_y")
    }
    
    # PRINTING OUTPUTS
    print(p3)
    
    # Optionally save the plot ---------------------
    
    if(!is.null(plot.file)){
      tmp <- utils.plot.save(p3,
                             dir=plot.dir,
                             file=plot.file,
                             verbose=verbose)
    }
    
    # # creating temp file names
    # if (plot.file) {
    #     temp_plot <- tempfile(pattern = "Plot_")
    #     match_call <-
    #         paste0(names(match.call()),
    #                "_",
    #                as.character(match.call()),
    #                collapse = "_")
    #     # saving to tempdir
    #     saveRDS(list(match_call, p3), file = temp_plot)
    #     if (verbose >= 2) {
    #         cat(report("  Saving the ggplot to session tempfile\n"))
    #     }
    # }
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    invisible(p3)
}
