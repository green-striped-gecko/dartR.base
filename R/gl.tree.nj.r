#' @name gl.tree.nj
#' @title Outputs a tree to summarize genetic similarity among populations (e.g. phenogram)
#' @family graphics

#' @description
#' This function is a wrapper for the nj function in package ape (and, for
#' method='upgma', the hclust function in stats) applied to Euclidean
#' distances calculated from the genlight object.

#' @details
#' A Euclidean distance matrix is calculated by default [dist.matrix = NULL]
#' from allele frequencies computed for each population at each locus,
#' excluding missing genotypes. Frequencies are scaled by ploidy, so SNP
#' (0/1/2) and Tag P/A (0/1) data both yield frequencies in the range 0 to 1.
#' Optionally the user can use as input for the tree any other distance
#' matrix using this parameter, see for example the function
#' \code{\link{gl.dist.pop}}.
#'
#' The tree is built with the neighbour-joining algorithm [method='nj'] or
#' by UPGMA clustering [method='upgma']; the legacy spelling 'ugpma' is
#' accepted as a synonym of 'upgma'. If an outgroup is specified, the tree
#' is rooted on that outgroup with the root resolved.
#'
#' The tree is computed and returned whether or not it is plotted; a failure
#' in plotting does not affect the returned tree.

#' @param x Name of the genlight object containing the SNP data [required].
#' @param dist.matrix Distance matrix [default NULL].
#' @param method Clustering method -- nj, neighbour-joining tree; upgma,
#' UPGMA tree [default 'nj'].
#' @param by.pop If TRUE, populations are the terminal taxa; if FALSE,
#' individuals are the terminal taxa [default TRUE].
#' @param as.pop Assign another ind.metric as the population for
#' the purposes of displaying more informative tip labels [default NULL].
#' @param type Type of dendrogram
#' "phylogram"|"cladogram"|"fan"|"unrooted"|"radial"|"tidy"
#'  [default "phylogram"].
#' @param outgroup Vector containing the population names that are the outgroups
#'  [default NULL].
#' @param labelsize Size of the labels as a proportion of the graphics default
#'  [default 0.7].
#' @param treefile Name of the file for the tree topology using Newick format
#' [default NULL].
#' @param plot.display If TRUE, the tree is plotted in the plot window
#' [default TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @return A tree of class phylo. If treefile is specified, the tree is also
#' written to that file in Newick format.
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # SNP data
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#'   gl.tree.nj(testset.gl,type='fan')
#' # Tag P/A data
#'   gl.tree.nj(testset.gs,type='fan')
#'   res <- gl.tree.nj(platypus.gl)
#'
#' @importFrom ape nj root plot.phylo write.tree
#' @importFrom graphics par
#' @importFrom stats hclust as.dist
#' @export

gl.tree.nj <- function(x,
                       dist.matrix = NULL,
                       method="nj",
                       by.pop=TRUE,
                       as.pop=NULL,
                       type = "phylogram",
                       outgroup = NULL,
                       labelsize = 0.7,
                       treefile = NULL,
                       plot.display = TRUE,
                       verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    if (verbose == 0) {
        plot.display <- FALSE
    }

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

    # FUNCTION SPECIFIC ERROR CHECKING
    method <- tolower(method)
    if (method == "ugpma") {
      # Legacy misspelling retained as a silent synonym for back compatibility
      method <- "upgma"
    }
    if (method != "nj" && method != "upgma") {
      if (verbose >= 1) {
        cat(warn("  Warning: method must be one of 'nj' or 'upgma'. Set to 'nj'. \n"))
      }
      method <- "nj"
    }

    if (!type %in% c("phylogram", "cladogram", "fan", "unrooted", "radial",
                     "tidy")) {
      stop(error("Fatal Error: type must be one of 'phylogram', 'cladogram',
                 'fan', 'unrooted', 'radial' or 'tidy'\n"))
    }

    # DO THE JOB

    if(by.pop==FALSE){
      pop(x) <- indNames(x)
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
        stop(error("Fatal Error: individual metric assigned to 'pop' does not exist. Check names(x@other$ind.metrics) and select again\n"))
      }
    }

    if(is.null(dist.matrix)){

      # Convert gl object to a matrix of allele frequencies, locus by population
      if (verbose >= 2) {
        cat(report(
          "  Converting to a matrix of frequencies, locus by populations\n"
        ))
      }
      # Frequencies exclude missing genotypes (na.rm), matching gl.dist.pop,
      # and are scaled by ploidy (2 for SNP, 1 for Tag P/A)
      pl <- ploidy(x)[1]
      t <- apply(as.matrix(x), 2, tapply, pop(x), function(e)
        mean(e, na.rm = TRUE) / pl)
      # Compute Euclidean distance
      if (verbose >= 2) {
        cat(report("  Computing Euclidean distances\n"))
      }
      d <- round(as.matrix(dist(t)), 4)
      d <- as.dist(d)
      # row.names(d) <- c(paste(row.names(d),' ')) row.names(d) <- substr(row.names(d),1,10)

    }else{
      # Coerce in case a plain matrix was supplied (hclust requires a dist object)
      d <- stats::as.dist(dist.matrix)
    }

    if(method=="upgma"){
      # Compute a UPGMA tree
      hc <- stats::hclust(d, method="average")
      tree <- as.phylo(hc)
    } else {
      # Compute an nj tree
      tree <- ape::nj(d)
    }

    if (!is.null(outgroup)) {
        # Root the tree on the outgroup, resolving the root node
        tree <- ape::root(tree, outgroup, resolve.root = TRUE)
    }

    # Output the tree file
    if (!is.null(treefile)) {
        if (verbose >= 2) {
            cat(report("  Writing the tree topology to", treefile, "\n"))
        }
        ape::write.tree(tree, file = treefile)
    }

    # Plot the tree; decoupled from the computation so a plotting failure
    # cannot lose the computed tree
    if (plot.display) {
        # Save the prior settings for mfrow, oma, mai and pty, and reassign
        op <-
            par(
                mfrow = c(1, 1),
                oma = c(1, 1, 1, 1),
                mai = c(0, 0, 0, 0),
                pty = "m"
            )
        on.exit(par(op))
        tryCatch(
            ape::plot.phylo(tree, type = type, cex = labelsize),
            error = function(e) {
                if (verbose >= 1) {
                    cat(warn("  Warning: tree computed but could not be plotted:",
                             conditionMessage(e), "\n"))
                }
            }
        )
    }

    # FLAG SCRIPT END

    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }

    return(tree)

}
