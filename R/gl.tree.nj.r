#' @name gl.tree.nj
#' @title Outputs a tree to summarize genetic similarity among populations (e.g. phenogram)
#' @family graphics

#' @description
#' This function is a wrapper for the nj function in package ape and hclust function in stats applied to Euclidean
#' distances calculated from the genlight object.

#' @param x Name of the genlight object containing the SNP data [required].
#' @param dist.matrix Distance matrix [default NULL].
#' @param method Clustering method -- nj, neighbor-joining tree; UGPMA, UGPMA tree [default 'nj'].
#' @param by.pop If TRUE, populations are the terminal taxa; if FALSE, individuals are the terminal taxa [default TRUE]
#' @param as.pop Assign another ind.metric as the population for
#' the purposes of displaying more informative tip labels [default NULL].
#' @param outgroup Vector containing the population names that are the outgroups
#'  [default NULL].
#' @param type Type of dendrogram "phylogram"|"cladogram"|"fan"|"unrooted"
#'  [default "phylogram"].
#' @param labelsize Size of the labels as a proportion of the graphics default
#'  [default 0.7].
#' @param treefile Name of the file for the tree topology using Newick format 
#' [default NULL].
#' @param verbose Verbosity: 0, silent, fatal errors only; 1, flag function
#' begin and end; 2, progress log; 3, progress and results summary; 5, full
#' report [default 2 or as specified using gl.set.verbosity].
#' 
#' @details
#' An euclidean distance matrix is calculated by default [dist.matrix = NULL]. 
#' Optionally the user can use as input for the tree any other distance matrix
#' using this parameter, see for example the function  \code{\link{gl.dist.pop}}.
#' 
#' @author Custodian: Arthur Georges (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#'  \donttest{
#' # SNP data
#'   gl.tree.nj(testset.gl,type='fan')
#' # Tag P/A data
#'   gl.tree.nj(testset.gs,type='fan')
#'   }
#'   res <- gl.tree.nj(platypus.gl)
#'   
#' @importFrom stringr str_pad
#' @importFrom ape nj root plot.phylo write.tree
#' @importFrom graphics hist par
#' @importFrom stats hclust
#' @export
#' @return A tree file of class phylo.

gl.tree.nj <- function(x,
                       dist.matrix = NULL,
                       method="nj",
                       by.pop=TRUE,
                       as.pop=NULL,
                       type = "phylogram",
                       outgroup = NULL,
                       labelsize = 0.7,
                       treefile = NULL,
                       verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    if (!is(x, "dartR")) {
      class(x) <- "dartR"  
      if (verbose>2) {
        cat(warn("Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"))
        cat(warn("                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR, 
                 please use function dartR2gl\n"))
      }
    }
    
    method <- tolower(method)
    if(method != "nj" && method != "ugpma"){
      cat(warn("  Warning: method must be one of nj or ugpma. Set to nj. \n"))
      method <- "nj"
    }
    
    # DO THE JOB
    
    if(by.pop==FALSE){
      popNames(x) <- IndNames(x)
      if (verbose >= 2) {
        cat(report("  Tree constructed for individuals\n"))
      }
    } else {
      if (verbose >= 2) {
        cat(report("  Tree constructed for populations\n"))
      }
    }
    
    # Assign the new population list if as.pop is specified -----------
    if (!is.null(as.pop)) {
      if (as.pop %in% names(x@other$ind.metrics)) {
        pop(x) <- unname(unlist(x@other$ind.metrics[as.pop]))
        if (verbose >= 2) {
          cat(report("  Assigning",as.pop,"as the tip labels\n"))
        }
      } else {
        stop(error("Fatal Error: individual metric assigned to 'pop' does not exist. Check names(gl@other$loc.metrics) and select again\n"))
      }
    }
    
    if(is.null(dist.matrix)){
      
      # Convert gl object to a matrix of allele frequencies, locus by population
      if (verbose >= 2) {
        cat(report(
          "  Converting to a matrix of frequencies, locus by populations\n"
        ))
      }
      t <- apply(as.matrix(x), 2, tapply, pop(x), function(e)
        mean(e) / 2)
      # Compute Euclidean distance
      if (verbose >= 2) {
        cat(report("  Computing Euclidean distances\n"))
      }
      d <- round(as.matrix(dist(t)), 4)
      d <- as.dist(d)
      # row.names(d) <- c(paste(row.names(d),' ')) row.names(d) <- substr(row.names(d),1,10)
      
    }else{
      d <- dist.matrix
    }
    
    if(method=="ugpma"){
      # Plot the distances as a UGPMA tree
      hc <- stats::hclust(d, method="average")
      tree <- as.phylo(hc)
    } else {
      # Plot the distances as an nj tree
      tree <- ape::nj(d)
    }
    
    if (!is.null(outgroup)) {
        # Function plot.phylo{ape} has the labels all of the same length outgroup <- stringr::str_pad(outgroup, nchar(tree$tip.label[1]),
        # side = c('right'), pad = ' ') # Truncate to 10 characters outgroup <- substr(outgroup,1,10) Root the tree
        tree <- ape::root(tree, outgroup)
        # Plot the tree Save the prior settings for mfrow, oma, mai and pty, and reassign
        op <-
            par(
                mfrow = c(1, 1),
                oma = c(1, 1, 1, 1),
                mai = c(0, 0, 0, 0),
                pty = "m"
            )
		on.exit(par(op))
        ape::plot.phylo(tree, type = type, cex = labelsize)
    } else {
        # Just plot the tree unrooted
        op <-
            par(
                mfrow = c(1, 1),
                oma = c(1, 1, 1, 1),
                mai = c(0, 0, 0, 0),
                pty = "m"
            )
        ape::plot.phylo(tree, type = type, cex = labelsize)
    }
    
    # Output the tree file
    if (!is.null(treefile)) {
        if (verbose >= 2) {
            cat(report("  Writing the tree topology to", treefile, "\n"))
        }
        ape::write.tree(tree, file = treefile)
    }
    
    
    # Reset the par options
	#now done by on exit
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(tree)
    
}
