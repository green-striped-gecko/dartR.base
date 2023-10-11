#' @name gl.smearplot
#' @title Smear plot of SNP or presence/absence (SilicoDArT) data
#' @family graphics

#' @description
#' Each locus is color coded for scores of 0, 1, 2 and NA for SNP data and 0, 1
#' and NA for presence/absence (SilicoDArT) data. Individual labels can be added.

#' Plot may become cluttered if ind.labels If there are too many individuals, 
#' it is best to use ind.labels = FALSE.
#' 
#' Works with both SNP data and P/A data (SilicoDArT)

#' @param x Name of the genlight object [required].
#' @param ind.labels If TRUE, individual IDs are shown [default FALSE].
# @param group.pop If ind.labels is TRUE, group by population [default TRUE].
#' @param label.size Size of the individual labels [default 10].
#' @param het.only If TRUE, show only the heterozygous state [default FALSE]
#' @param plot.display If TRUE, the plot is displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors List of four color names for the column fill for homozygous reference,
#' heterozygous, homozygous alternate, and missing value (NA) [default c("#0000FF","#00FFFF","#FF0000","#e0e0e0")].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]#' 
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param legend Position of the legend: “left”, “top”, “right”, “bottom” or
#'  'none' [default = 'bottom'].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity]
#' 
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' gl.smearplot(testset.gl,ind.labels=FALSE)
#' gl.smearplot(testset.gs,ind.labels=FALSE)
#' gl.smearplot(testset.gl[1:10,],ind.labels=TRUE)
#' gl.smearplot(testset.gs[1:10,],ind.labels=TRUE)

#' @export
#' @return Returns the ggplot object
#' 
# TEST
# ddd <- matrix(data=0,nrow=10,ncol=10)
# ddd[8,10] <- NA
# ddd[9,10] <- 2
# ddd[10,10] <- 2
# ddd
# ddd <- as.genlight(ddd)
# ploidy(ddd) <- 2
# ddd <- gl.compliance.check(ddd)
# gl.smearplot(ddd)

gl.smearplot <- function(x,
                         plot.display=TRUE,
                         ind.labels = FALSE,
                         label.size = 10,
                         #group.pop = FALSE, 
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
  if(verbose==0){plot.display <- FALSE}
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir,verbose=0)
  
  # SET COLOURS
  if(is.null(plot.colors)){
    plot.colors <- c("#0000FF","#00FFFF","#FF0000","#e0e0e0")
  } else {
    if(length(plot.colors)>4){
      if(verbose >= 2)cat(warn("  Specified plot colours exceed 4, first 4 only are used\n"))
      plot.colors <- plot.colors[1:4]
    }
  }
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose) 
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)
  
  if (het.only) {
    plot.colors <- c("#d3d3d3","#00FFFF","#d3d3d3","#e0e0e0")
  }
  
  # SET IND LABELS
  
  if(ind.labels == TRUE){
    individuals <- indNames(x)
  } else {
    individuals <- seq(1:nInd(x))
  }
  
  # DO THE JOB
  
  # pull the data from the genlight object, and place in a dataframe
  df.matrix <- as.data.frame(as.matrix(x))
  colnames(df.matrix) <- 1:nLoc(x)
  df.matrix$id <- individuals
  df.matrix$pop <- pop(x)
  
  # convert the data to long form
  df.listing <- reshape2::melt(df.matrix, id.vars = c("pop", "id"))
  df.listing$value <- as.character(df.listing$value)
  df.listing$value <- ifelse(df.listing$value=="NA", NA, df.listing$value)
  colnames(df.listing) <- c("pop", "id", "locus", "genotype")
  df.listing$id <- as.factor(df.listing$id)
  
  # The locus names are 1 to nLoc(x)
  
  loc_labels <- pretty(1:nLoc(x), 5)
  id_labels <- pretty(1:nInd(x), 5)
  
  locus <- id <- genotype <- NA
  
  # Assign colours and labels for genotypic data
  labels_genotype <- as.character(unique(df.listing$genotype))
  labels_genotype <- labels_genotype[!is.na(labels_genotype)]
  labels_genotype <- labels_genotype[order(labels_genotype)]
  #labels_genotype <- c("0","1","2")
  plot.colors.hold <- plot.colors
  tmp <- NULL
  if(length(labels_genotype) < 3){
    if("0" %in% labels_genotype){
      tmp[1] <- plot.colors[1]
    }
    if ("1" %in% labels_genotype){
      if(is.null(tmp)){
        tmp <- plot.colors[2]
      } else {
        tmp <- c(tmp,plot.colors[2])
      }
    }
    if ("2" %in% labels_genotype){
      if(is.null(tmp)){
        tmp <- plot.colors[3]
      } else {
        tmp <- c(tmp,plot.colors[3])
      }
    }
    tmp <- c(tmp,plot.colors[4])
    
    plot.colors <- tmp
  }
  n.colors <- length(plot.colors)
  
  labels_genotype[which(is.na(labels_genotype))] <- "Missing data"
  labels_genotype[labels_genotype=="0"] <- "Homozygote reference"
  labels_genotype[labels_genotype=="1"] <- "Heterozygote"
  labels_genotype[labels_genotype=="2"] <- "Homozygote alternate"
  
  
  if (datatype == "SNP") {
    if(ind.labels==TRUE){
      p3 <-
        ggplot(df.listing, aes(
          x = locus,
          y = id,
          fill = genotype
        )) + geom_raster() + 
        scale_fill_discrete(
          type = plot.colors,
          na.value = plot.colors[n.colors],
          name = "Genotype",
          labels = labels_genotype
          # ) + theme_dartR() + theme(
        ) + theme(
          legend.position = legend,
          #axis.text.y = element_text(size = label.size)
          axis.text.y = element_text(size = label.size)
        ) +
        scale_x_discrete(
          breaks = loc_labels,
          labels = as.character(loc_labels),
          name = "Loci",
          position="bottom"
        ) +
        ylab("Individuals")
    } else {
      p3 <-
        ggplot(df.listing, aes(
          x = locus,
          y = id,
          fill = genotype
        )) + geom_raster() + 
        scale_fill_discrete(
          type = plot.colors,
          na.value = plot.colors[n.colors],
          name = "Genotype",
          labels = labels_genotype
          # ) + theme_dartR() + theme(
        ) + theme(
          legend.position = legend,
          #axis.text.y = element_text(size = label.size)
          axis.text.y = element_text(size = label.size)
        ) +
        scale_x_discrete(
          breaks = loc_labels,
          labels = as.character(loc_labels),
          name = "Loci",
          position="bottom"
        ) +
        scale_y_discrete(
          breaks = id_labels,
          labels = as.character(id_labels),
          name = "Individuals",
          position="left"
        )
      #ylab("Individuals")
    }
  }
  
  # Assign labels for presence absence data
  #labels_silicodart <- as.character(unique(df.listing$genotype))
  labels_silicodart <- c("0","1")
  labels_silicodart[which(is.na(labels_silicodart))] <- "Missing data"
  labels_silicodart["0"] <- "Absence"
  labels_silicodart["1"] <- "Presence"
  
  plot.colors <- plot.colors.hold
  
  if (datatype == "SilicoDArT") {
    if(het.only){
      cat(warn("The het only option is applicable to SNP data only. Set to FALSE\n"))
      het.only <- FALSE
    }
    if(ind.labels==TRUE){
      p3 <-
        ggplot(df.listing, aes(
          x = locus,
          y = id,
          fill = genotype
        )) + geom_raster() + scale_fill_discrete(
          type = plot.colors[c(1,3)],
          na.value = plot.colors[4],
          name = "Sequence Tag",
          labels = labels_silicodart
        ) + theme_dartR() + theme(
          legend.position = legend,
          axis.text.y = element_text(size = label.size)
        ) +
        scale_x_discrete(
          breaks = loc_labels,
          labels = as.character(loc_labels),
          name = "Loci"
        ) +
        ylab("Individuals")
    } else {
      p3 <-
        ggplot(df.listing, aes(
          x = locus,
          y = id,
          fill = genotype
        )) + geom_raster() + scale_fill_discrete(
          type = plot.colors[c(1,3)],
          na.value = plot.colors[4],
          name = "Sequence Tag",
          labels = labels_silicodart
        ) + theme_dartR() + theme(
          legend.position = legend,
          axis.text.y = element_text(size = label.size)
        ) +
        scale_x_discrete(
          breaks = loc_labels,
          labels = as.character(loc_labels),
          name = "Loci"
        ) +
        scale_y_discrete(
          breaks = id_labels,
          labels = as.character(id_labels),
          name = "Individuals",
          position="left"
        )
      #ylab("Individuals")
    }
  }
  
  
  # if (ind.labels==TRUE & group.pop == TRUE) {
  #     p3 <- p3 + facet_wrap(~ pop,
  #                           ncol = 1,
  #                           dir = "v",
  #                           scales = "free_y")
  # }
  
  # PRINTING OUTPUTS
  print(p3)
  
  # Optionally save the plot ---------------------
  
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
  
  invisible(p3)
}
