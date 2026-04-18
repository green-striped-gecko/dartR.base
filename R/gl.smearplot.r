#' @name gl.smearplot
#' @title A plot of individuals against loci, showing the character state [0, 1, 2, NA]
#' @family graphics

#' @description
#' Each locus is color coded for scores of 0, 1, 2 and NA for SNP data and 0, 1
#' and NA for presence/absence (SilicoDArT) data. Individual labels can be added.
#'
#' Plot may become cluttered if individual labels are shown and there are many
#' individuals. In such cases, it is best to use `ind.labels = FALSE`.
#'
#' Works with both SNP data and presence/absence data (SilicoDArT).
#'
#' @param x Name of the genlight object [required].
#' @param ind.labels If TRUE, individual IDs are shown [default FALSE].
#' @param ind.order Optional character vector of length nInd(x) that determines the display order of individuals
#'   on the y-axis. If `ind.order` is used, the accompanying tree will not be
#'   displayed.
#' @param loc.order Optional numeric vector of length `nLoc(x)` giving a locus
#'   metric used to order loci on the x-axis (for example
#'   `x@other$locmetrics$avgPIC`). Loci are plotted in increasing order of this
#'   vector; use `-locmetric` for decreasing order.
#' @param group.pop If TRUE, facet the plot by population [default FALSE].
#' @param label.size Size of the individual labels [default 10].
#' @param het.only If TRUE, show only the heterozygous state [default FALSE].
#' @param plot.display If TRUE, the plot is displayed in the plot window
#'   [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#'   [default NULL].
#' @param plot.colors List of four colour values for homozygous reference,
#'   heterozygous, homozygous alternate, and missing value (NA)
#'   [default c("#0000FF", "#00FFFF", "#FF0000", "#e0e0e0")].
#' @param plot.dir Directory in which to save the plot files
#'   [default as specified by the global working directory or `tempdir()`].
#' @param plot.file Base name for the plot file to save (exclude extension)
#'   [default NULL].
#' @param legend Position of the legend: `"left"`, `"top"`, `"right"`,
#'   `"bottom"` or `"none"` [default `"bottom"`].
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
#' @return Returns the ggplot object invisibly

gl.smearplot <- function(x,
                        plot.display=TRUE,
                        ind.labels = FALSE,
                        ind.order = NULL,
                        loc.order = NULL,
                        label.size = 10,
                        group.pop = FALSE, 
                        plot.theme = NULL,
                        plot.colors = NULL,
                        plot.file=NULL,
                        plot.dir=NULL,
                        het.only=FALSE,
                        legend = "bottom",
                        verbose = NULL) {
    
  # CHECK IF PACKAGES ARE INSTALLED
  pkg <- "reshape2"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Package '", pkg,
        "' is required for gl.smearplot(). Please install it."
      ),
      call. = FALSE
    )
  }
    
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    if(verbose==0){plot.display <- FALSE}
    
    # SET WORKING DIRECTORY
    plot.dir <- gl.check.wd(plot.dir,verbose=0)
    
    # SET COLOURS
    default.plot.colors <- c("#0000FF", "#00FFFF", "#FF0000", "#e0e0e0")
    
    if (is.null(plot.colors)) {
      plot.colors <- default.plot.colors
    } else {
      if (length(plot.colors) < 4) {
        if (verbose >= 2) {
          cat(warn(
            "  Specified plot colours fewer than 4; missing values replaced with defaults\n"
          ))
        }
        plot.colors <- c(plot.colors, default.plot.colors[(length(plot.colors) + 1):4])
      } else if (length(plot.colors) > 4) {
        if (verbose >= 2) {
          cat(warn(
            "  Specified plot colours exceed 4, first 4 only are used\n"
          ))
        }
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
    
    # CHECK AND SET INDIVIDUAL ORDER -----------------------------------
    
    ind.names <- indNames(x)
    use.ind.order <- !is.null(ind.order)
    
    if (use.ind.order) {
      if (is.null(ind.names)) {
        stop(
          "ind.order cannot be used because indNames(x) is NULL.",
          call. = FALSE
        )
      }
      
      if (!is.character(ind.order)) {
        stop(
          "ind.order must be a character vector giving a permutation of indNames(x).",
          call. = FALSE
        )
      }
      
      if (length(ind.order) != nInd(x)) {
        stop(
          paste0(
            "ind.order must have length ", nInd(x),
            " to match the number of individuals in x."
          ),
          call. = FALSE
        )
      }
      
      if (anyDuplicated(ind.order)) {
        stop(
          "ind.order contains duplicates; it must be a permutation of indNames(x).",
          call. = FALSE
        )
      }
      
      if (!setequal(ind.order, ind.names)) {
        missing.names <- setdiff(ind.names, ind.order)
        extra.names <- setdiff(ind.order, ind.names)
        
        msg <- "ind.order must be a permutation of indNames(x)."
        if (length(missing.names) > 0) {
          msg <- paste0(
            msg, "\nMissing from ind.order: ",
            paste(missing.names, collapse = ", ")
          )
        }
        if (length(extra.names) > 0) {
          msg <- paste0(
            msg, "\nNot found in indNames(x): ",
            paste(extra.names, collapse = ", ")
          )
        }
        
        stop(msg, call. = FALSE)
      }
      
      if (verbose >= 1) {
        cat(report("Using user-specified individual order from ind.order\n"))
      }
      
      # Replace 'plot.tree' below with the name of your actual tree flag
      if (exists("plot.tree", inherits = FALSE) && isTRUE(plot.tree)) {
        plot.tree <- FALSE
        if (verbose >= 1) {
          cat(warn(
            "  ind.order supplied: accompanying tree will not be displayed\n"
          ))
        }
      } else {
        if (verbose >= 1) {
          cat(report(
            "ind.order supplied: any accompanying tree should be suppressed\n"
          ))
        }
      }
      
      row.order <- match(ind.order, ind.names)
      
    } else {
      row.order <- seq_len(nInd(x))
    }
    
    # SET IND LABELS
    
    if (ind.labels == TRUE) {
      individuals <- ind.names[row.order]
    } else {
      individuals <- seq_len(nInd(x))
    }
    
    # CHECK AND SET LOCUS ORDER ----------------------------------------
    
    if (!is.null(loc.order)) {
      if (!is.numeric(loc.order)) {
        stop(
          "loc.order must be a numeric vector of length nLoc(x).",
          call. = FALSE
        )
      }
      
      if (length(loc.order) != nLoc(x)) {
        stop(
          paste0(
            "loc.order must have length ", nLoc(x),
            " to match the number of loci in x."
          ),
          call. = FALSE
        )
      }
      
      if (all(is.na(loc.order))) {
        stop(
          "loc.order contains only NA values; loci cannot be ordered.",
          call. = FALSE
        )
      }
      
      col.order <- order(loc.order, na.last = TRUE)
      
      if (verbose >= 1) {
        cat(report(
          "Ordering loci on the x-axis by loc.order (ascending; NA last)\n"
        ))
      }
      if (anyNA(loc.order) && verbose >= 2) {
        cat(warn(
          "  loc.order contains NA values; corresponding loci placed at the right-hand end\n"
        ))
      }
    } else {
      col.order <- seq_len(nLoc(x))
    }
    
    display.loc.ids <- as.character(seq_len(nLoc(x)))
    
    # DO THE JOB
    
    # pull the data from the genlight object, reorder if requested,
    # and place in a dataframe
    df.matrix <- as.data.frame(as.matrix(x))
    df.matrix <- df.matrix[row.order, col.order, drop = FALSE]
    colnames(df.matrix) <- display.loc.ids
    df.matrix$id <- individuals
    df.matrix$pop <- pop(x)[row.order]
    
    # convert the data to long form
    df.listing <- reshape2::melt(df.matrix, id.vars = c("pop", "id"))
    df.listing$value <- as.character(df.listing$value)
    df.listing$value <- ifelse(df.listing$value=="NA", NA, df.listing$value)
    colnames(df.listing) <- c("pop", "id", "locus", "genotype")
    df.listing$id <- factor(df.listing$id,levels = unique(df.matrix$id))
    df.listing$locus <- factor(df.listing$locus, levels = display.loc.ids)
    
    # Number loci on the x-axis from 1 to nLoc(x), regardless of ordering
    
    if (nLoc(x) <= 20) {
      loc_breaks <- display.loc.ids
      loc_labels <- display.loc.ids
    } else {
      loc_idx <- pretty(seq_len(nLoc(x)), 5)
      loc_idx <- unique(as.integer(loc_idx))
      loc_idx <- loc_idx[loc_idx >= 1 & loc_idx <= nLoc(x)]
      loc_breaks <- as.character(loc_idx)
      loc_labels <- as.character(loc_idx)
    }
    
    id_labels <- pretty(seq_len(nInd(x)), 5)
    
    locus <- id <- genotype <- NA
    
    labels_genotype <- as.character(unique(df.listing$genotype))
    has_na <- any(is.na(df.listing$genotype))
    
    labels_genotype_no_na <- labels_genotype[!is.na(labels_genotype)]
    labels_genotype_no_na <- labels_genotype_no_na[order(labels_genotype_no_na)]
    
    plot.colors.hold <- plot.colors
    tmp <- NULL
    if(length(labels_genotype_no_na) < 3){
      if("0" %in% labels_genotype_no_na){
        tmp[1] <- plot.colors[1]
      }
      if("1" %in% labels_genotype_no_na){
        if(is.null(tmp)){
          tmp <- plot.colors[2]
        } else {
          tmp <- c(tmp, plot.colors[2])
        }
      }
      if("2" %in% labels_genotype_no_na){
        if(is.null(tmp)){
          tmp <- plot.colors[3]
        } else {
          tmp <- c(tmp, plot.colors[3])
        }
      }
      if(has_na){
        tmp <- c(tmp, plot.colors[4])
      }
      plot.colors <- tmp
    }
    n.colors <- length(plot.colors)
    
    labels_genotype <- labels_genotype_no_na
    labels_genotype[labels_genotype == "0"] <- "Homozygote reference"
    labels_genotype[labels_genotype == "1"] <- "Heterozygote"
    labels_genotype[labels_genotype == "2"] <- "Homozygote alternate"
    if(has_na){
      labels_genotype <- c(labels_genotype, "Missing data")
    }
    

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
          ) + 
          scale_x_discrete(
            breaks = loc_breaks,
            labels = loc_labels,
            name = "Loci",
            position="bottom"
          ) +
        ylab("Individuals")
        
        if(!is.null(plot.theme)){
          p3 <- p3 + plot.theme
        }
        
        p3 <- p3  + theme(
          legend.position = legend,
          #axis.text.y = element_text(size = label.size)
          axis.text.y = element_text(size = label.size)
        ) 
        
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
            ) + 
            scale_x_discrete(
              breaks = loc_breaks,
              labels = loc_labels,
                name = "Loci",
                position="bottom"
            ) +
            scale_y_discrete(
              breaks = as.character(id_labels),
              labels = as.character(id_labels),
                name = "Individuals",
                position="left"
          )
        #ylab("Individuals")
        
        if(!is.null(plot.theme)){
          p3 <- p3 + plot.theme
        }
        
        p3 <- p3  + theme(
          legend.position = legend,
          #axis.text.y = element_text(size = label.size)
          axis.text.y = element_text(size = label.size)
        ) 
      }
    }
    
    # Assign labels for presence absence data
    has_na_silicodart <- any(is.na(df.listing$genotype))
    
    labels_silicodart <- c("Absence", "Presence")
    if (has_na_silicodart) {
      labels_silicodart <- c(labels_silicodart, "Missing data")
    }
    
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
            na.translate = TRUE,
            name = "Sequence Tag",
            labels = labels_silicodart
          ) + 
          scale_x_discrete(
            breaks = loc_breaks,
            labels = loc_labels,
            name = "Loci"
          ) +
          ylab("Individuals")
        
        if(!is.null(plot.theme)){
          p3 <- p3 + plot.theme
        }
        
        p3 <- p3  + theme(
          legend.position = legend,
          axis.text.y = element_text(size = label.size)
        ) 
        
      } else {
        p3 <-
          ggplot(df.listing, aes(
            x = locus,
            y = id,
            fill = genotype
          )) + geom_raster() + scale_fill_discrete(
            type = plot.colors[c(1,3)],
            na.value = plot.colors[4],
            na.translate = TRUE,
            name = "Sequence Tag",
            labels = labels_silicodart
          ) + 
          scale_x_discrete(
            breaks = loc_breaks,
            labels = loc_labels,
            name = "Loci"
          ) +
          scale_y_discrete(
            breaks = as.character(id_labels),
            labels = as.character(id_labels),
            name = "Individuals",
            position="left"
          )
        
        if(!is.null(plot.theme)){
          p3 <- p3 + plot.theme
        }
        
        p3 <- p3  + theme(
          legend.position = legend,
          axis.text.y = element_text(size = label.size)
        ) 
        
      }
    }
    
    if (group.pop == TRUE) {
        p3 <- p3 + facet_wrap(~ pop,
                              ncol = 1,
                              dir = "v",
                              scales = "free_y")
    }
    
    # PRINTING OUTPUTS
    if (plot.display) print(p3)
    
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
