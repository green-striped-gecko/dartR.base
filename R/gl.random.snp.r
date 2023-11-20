#' @name gl.randomize.snps
#' @title Randomly changes the allocation of 0's and 2's in a genlight object
#' @description
#' This function samples randomly half of the SNPs and re-codes, in the sampled
#' SNP's, 0's by 2's.

#' @param x Name of the genlight object containing the SNP data [required].
#' @param plot.display If TRUE, resultant plots are displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors List of two color names for the borders and fill of the
#'  plots [default c("#2171B5","#6BAED6")].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report [default NULL,
#' unless specified using gl.set.verbosity].

#' @details
#' DArT calls the most common allele as the reference allele. In a genlight
#' object, homozygous for the reference allele are coded with a '0' and
#' homozygous for the alternative allele are coded with a '2'. This causes some
#' distortions in visuals from time to time.

#' If plot.display = TRUE, two smear plots (pre-randomisation and
#' post-randomisation) are presented using a random subset of individuals (10)
#' and loci (100) to provide an overview of the changes.

#' Resultant ggplots are saved to the session's temporary directory.

#' @return Returns a genlight object with half of the loci re-coded.
#' @author Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' require("dartR.data")
#' res <- gl.randomize.snps(platypus.gl[1:5,1:5],verbose = 5)

#' @export

gl.randomize.snps <- function(x,
                              plot.display=TRUE,
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
    plot.colors <- c("#0000FF","#00FFFF","#FF0000","#e0e0e0")
  } else {
    if(length(plot.colors) > 4){
      if(verbose >= 2){cat(warn("  More than 2 colors specified, only the first 4 are used\n"))}
      plot.colors <- plot.colors[1:4]
    }
  }
  
  # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # DO THE JOB
    
    hold <- x
    
    snp_matrix_temp <- as.matrix(x)
    snp_matrix_temp_0 <- snp_matrix_temp == 0
    snp_matrix_temp_2 <- snp_matrix_temp == 2
    
    snp_matrix_temp[snp_matrix_temp_0 == TRUE] <- 2
    snp_matrix_temp[snp_matrix_temp_2 == TRUE] <- 0
    
    random_snps <- sample(1:nLoc(x), nLoc(x) / 2)
    
    snp_matrix <- as.matrix(x)
    snp_matrix[, random_snps] <- snp_matrix_temp[, random_snps]
    
    x@gen <-
        lapply(1:nrow(snp_matrix), function(i)
            new("SNPbin", as.integer(snp_matrix[i, ])))
    
    random_snps <- random_snps[order(random_snps)]
    
    if (verbose == 5) {
        cat(report(paste(
            "The loci that were changed are:",
            paste(random_snps, collapse = ", "),
            "\n"
        )))
    }
    
    if (plot.display) {
        # subsetting objects to provide an overview of the changes
        if (nInd(x) > 10) {
            ind_to_plot <- sample(1:nInd(x), 10)
            x_plot <- x[ind_to_plot, ]
            hold_plot <- hold[ind_to_plot, ]
        } else {
            x_plot <- x
            hold_plot <- hold
        }
        if (nLoc(x_plot) > 100) {
            loc_to_plot <- sample(1:nLoc(x_plot), 100)
            x_plot <- x_plot[, loc_to_plot]
            hold_plot <- hold_plot[, loc_to_plot]
        }
        
        # plot before randomisation
        p1 <-
            gl.smearplot(hold_plot, legend = "none", verbose = 0)
        p1 <-
            p1 + ggtitle("Pre-randomisation") + theme(
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
            )
        
        # plot after randomisation
        p2 <- gl.smearplot(x_plot, verbose = 0)
        p2 <- p2 + ggtitle("Post-randomisation")
    }
    
    # PRINTING OUTPUTS
    
      # using package patchwork
        p3 <- p1 / p2
        if (plot.display) {print(p3)}
    
        # Optionally save the plot ---------------------
        
        if(!is.null(plot.file)){
          tmp <- utils.plot.save(p3,
                                 dir=plot.dir,
                                 file=plot.file,
                                 verbose=verbose)
        }
    
    # ADD TO HISTORY
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    invisible(x)
    
}
