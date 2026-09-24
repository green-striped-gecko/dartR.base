#' @name gl.randomize.snps
#' @title Randomly swaps the homozygote coding (0 and 2) in half of the loci
#' @description
#' This function samples half of the loci at random and, in the sampled loci,
#' swaps the coding of the two homozygotes: 0 becomes 2 and 2 becomes 0.
#' Heterozygotes (1) and missing values are unchanged. The allele labels of the
#' sampled loci are reversed to match, so the object still describes the same
#' genotypes; only which allele is counted changes.

#' @param x Name of the genlight object containing the SNP data [required].
#' @param plot.display If TRUE, resultant plots are displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors Vector of four color names passed to
#' \code{\link{gl.smearplot}} for homozygote reference, heterozygote,
#' homozygote alternative and missing data
#' [default c("#0000FF","#00FFFF","#FF0000","#e0e0e0")].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL]
#' @param plot.dir Directory in which to save the plot file [default as
#' specified by the global working directory or tempdir()]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report [default NULL,
#' unless specified using gl.set.verbosity].

#' @details
#' DArT calls the most common allele as the reference allele. In a genlight
#' object, homozygous for the reference allele are coded with a '0' and
#' homozygous for the alternative allele are coded with a '2'. This causes some
#' distortions in visuals from time to time.

#' Only SNP data are accepted. FBM-backed objects are recoded in a copy of
#' their backing file, so the input object is not modified.

#' If plot.display = TRUE, two smear plots (pre-randomisation and
#' post-randomisation) are presented using a random subset of individuals (10)
#' and loci (100) to provide an overview of the changes. The combined plot is
#' saved only when plot.file is specified.

#' @return Returns a genlight object with half of the loci re-coded and their
#' allele labels (loc.all) reversed. Locus metric flags are reset, as the
#' allele-frequency metrics no longer match the recoded genotypes.
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' gl <- gl.filter.monomorphs(testset.gl)
#' res <- gl.randomize.snps(gl,verbose = 5)

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
      if(verbose >= 2){cat(warn("  More than 4 colors specified, only the first 4 are used\n"))}
      plot.colors <- plot.colors[1:4]
    }
  }

  # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)

    # CHECK DATATYPE
    # Presence/absence data have no second homozygote; swapping would write
    # invalid 2 codes into a ploidy-1 object
    datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

    # DO THE JOB

    # Keep the unrecoded object only for the before/after plot
    if (plot.display) {
      hold <- x
    }

    random_snps <- sort(sample(seq_len(nLoc(x)), floor(nLoc(x) / 2)))

    # 2 - g swaps the homozygotes (0 <-> 2) and leaves heterozygotes (1) and
    # NA unchanged
    fbm <- .fbm_or_null(x)
    if (is.null(fbm)) {
      snp_matrix <- as.matrix(x)
      snp_matrix[, random_snps] <- 2 - snp_matrix[, random_snps]
      x@gen <- matrix2gen(snp_matrix)
    } else {
      # FBM storage is shared by reference: writing into x@fbm would also
      # recode the caller's object, so the recoding goes into a copy. Genotype
      # reads come from the FBM, not @gen, so @gen is left as it is. NA is
      # stored as code 3 (CODE_012 mapping, as in gl.gen2fbm). Columns are
      # recoded in blocks to avoid densifying the whole matrix.
      x@fbm <- bigstatsr::big_copy(fbm, backingfile = tempfile("geno_"))
      blocks <- split(random_snps, ceiling(seq_along(random_snps) / 1000))
      for (cols in blocks) {
        block <- 2 - x@fbm[, cols, drop = FALSE]
        block[is.na(block)] <- 3
        x@fbm[, cols] <- block
      }
    }

    # Swapping the homozygote codes switches which allele is counted, so the
    # allele labels are reversed to keep describing the same genotypes
    if (length(x@loc.all) == nLoc(x)) {
      x@loc.all[random_snps] <-
        vapply(strsplit(x@loc.all[random_snps], "/", fixed = TRUE),
               function(alleles) paste(rev(alleles), collapse = "/"),
               character(1))
    }

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
            gl.smearplot(hold_plot, legend = "none", plot.theme = plot.theme,
                         plot.colors = plot.colors, verbose = 0)
        p1 <-
            p1 + ggtitle("Pre-randomisation") + theme(
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
            )

        # plot after randomisation
        p2 <- gl.smearplot(x_plot, plot.theme = plot.theme,
                           plot.colors = plot.colors, verbose = 0)
        p2 <- p2 + ggtitle("Post-randomisation")
    }

    # PRINTING OUTPUTS
    # p1 and p2 are built only inside the plot.display block above, so assemble,
    # print and (optionally) save the combined plot only when plot.display =
    # TRUE. Otherwise `p3 <- p1 / p2` errors with "object 'p1' not found" -- for
    # plot.display = FALSE and for verbose = 0, which forces plot.display FALSE.
    if (plot.display) {
      # using package patchwork
        p3 <- p1 / p2
        print(p3)

        # Optionally save the plot ---------------------

        if(!is.null(plot.file)){
          tmp <- utils.plot.save(p3,
                                 dir=plot.dir,
                                 file=plot.file,
                                 verbose=verbose)
        }
    }

    # RESET FLAGS
    # Half the loci had their 0/2 homozygote coding swapped, so the per-locus
    # allele-frequency metrics (OneRatioRef/OneRatioSnp, FreqHomRef/FreqHomSnp,
    # PICRef/PICSnp, maf, ...) no longer match the recoded genotypes. Mark them
    # as no longer current so downstream gl.recalc.metrics / compliance checks
    # know to recompute them.
    x <- utils.reset.flags(x, verbose = 0)

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
