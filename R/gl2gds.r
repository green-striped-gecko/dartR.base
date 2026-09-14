#' @name gl2gds
#' @title Converts a genlight object into gds format
#' @family linker

#' @description
#' Package SNPRelate relies on a bit-level representation of a SNP dataset that
#' competes with \{adegenet\} genlight objects and associated files. This
#' function converts a genlight object to a gds format file.
#' 
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file (including extension)
#' [default 'gl_gds.gds'].
#' @param outpath Path where to save the output file [default global working 
#' directory or if not specified, tempdir()].
#' @param snp.pos Field name from the slot loc.metrics where the SNP position is
#' stored [default '0'].
#' @param snp.chr Field name from the slot loc.metrics where the chromosome of
#' each is stored [default '0'].
#' @param chr.format Whether chromosome information is stored as 'numeric' or as
#' 'character', see details [default 'character'].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#' 
#' @details
#' This function orders the SNPS by chromosome and by position before converting
#' to SNPRelate format, as required by this package.

#' The chromosome of each SNP can be a character or numeric, as described in the
#' vignette of SNPRelate:
#' 'snp.chromosome, an integer or character mapping for each chromosome.
#' Integer: numeric values 1-26, mapped in order from 1-22, 23=X, 24=XY
#' (the pseudoautosomal region), 25=Y, 26=M (the mitochondrial probes), and 0
#' for probes with unknown positions; it does not allow NA. Character: “X”,
#'  “XY”, “Y” and “M” can be used here, and a blank string indicating unknown
#'  position.'

#' When using some functions from package SNPRelate with datasets other than
#' humans it might be necessary to use the option autosome.only=FALSE to avoid
#' detecting chromosome coding. So, it is important to read the documentation of
#' the function before using it.

#' The chromosome information for unmapped SNPS is coded as 0, as required by
#' SNPRelate.

#' Remember to close the GDS file before working in a different GDS object with
#' the function \link[SNPRelate]{snpgdsClose} (package SNPRelate).
#' 
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' \donttest{
#' require("dartR.data")
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' gl2gds(platypus.gl,snp.pos='ChromPos_Platypus_Chrom_NCBIv1',
#' snp.chr = 'Chrom_Platypus_Chrom_NCBIv1', outpath=tempdir())
#' }
#' 
#' @importFrom SNPRelate snpgdsCreateGeno snpgdsOpen snpgdsSummary snpgdsClose
#' @export
#' @return  returns no value (i.e. NULL)

gl2gds <- function(x,
                   outfile = "gl_gds.gds",
                   outpath = NULL,
                   snp.pos = "0",
                   snp.chr = "0",
                   chr.format = "character",
                   verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # SET WORKING DIRECTORY
    outpath <- gl.check.wd(outpath,verbose=0)
     outfilespec <- file.path(outpath, outfile)
     
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

    # DO THE JOB

    # ordering loc.metrics by chromosome and snp position
    snp_order_temp <- x$other$loc.metrics
    snp_order_temp$snp_id <- locNames(x)

    if (snp.chr == 0) {
        snp_order_temp$chrom <- 0
    } else {
        if (chr.format == "numeric") {
            snp_order_temp$chrom <-
                as.numeric(unname(unlist(snp_order_temp[snp.chr])))
        }
        if (chr.format == "character") {
            snp_order_temp$chrom <-
                as.character(unname(unlist(snp_order_temp[snp.chr])))
        }
    }
    
    if (snp.pos == 0) {
        snp_order_temp$snp.pos <- 0
    } else {
        snp_order_temp$snp.pos <-
            as.numeric(unname(unlist(snp_order_temp[snp.pos])))
    }
    
    # Convert any NA values to 0 (genlight objects have NA for missing; 
    #SNPRelate has 0 in this instance)
    snp_order_temp[is.na(snp_order_temp$snp.pos), "snp.pos"] <-
        0
    # Convert any NA values to 0 (genlight objects have NA for missing; 
    #SNPRelate has 0 in this instance)
    snp_order_temp[snp_order_temp$snp.pos == 0, "chrom"] <- 0

    # one sort permutation, applied to EVERY per-locus field, so genotype
    # rows, snp.id, snp.allele, positions and chromosomes stay aligned
    # record by record
    ord <- order(snp_order_temp$chrom, snp_order_temp$snp.pos)

    # ordering snp matrix; SNPRelate defines the stored genotype as the
    # count of the FIRST allele of snp.allele, while the genlight dosage
    # counts the second, so the dosage is inverted (2 - dosage) with
    # missing coded 3
    genmat_temp <- t(as.matrix(x))
    genmat_temp <- genmat_temp[ord, , drop = FALSE]
    genmat_temp <- 2 - genmat_temp
    genmat_temp[is.na(genmat_temp)] <- 3

    snp.id_temp <- locNames(x)[ord]

    snp.allele_temp <- x@loc.all[ord]

    sample.id_temp <- indNames(x)
    sample.id_temp <-
        gsub(" ", replacement = "_", sample.id_temp)

    geno_list <-
        list(
            sample.id = sample.id_temp,
            snp.id = snp.id_temp,
            snp.position = snp_order_temp$snp.pos[ord],
            snp.chromosome = snp_order_temp$chrom[ord],
            snp.allele = snp.allele_temp,
            genotype = genmat_temp
        )
    
    # Create the gds file
    if (verbose >= 2) {
        cat(report("  Converting SNP data to gds format\n"))
    }
    
    # create a gds file
    with(
        geno_list,
        SNPRelate::snpgdsCreateGeno(
            gds.fn = outfilespec,
            genmat = genotype,
            sample.id = sample.id,
            snp.id = snp.id,
            snp.chromosome = snp.chromosome,
            snp.position = snp.position,
            snp.allele = snp.allele,
            snpfirstdim = TRUE
        )
    )
    
    # Open the GDS file, which will print out a summary of contents
    if (verbose >= 2) {
        cat(report("  Writing data to file", outfilespec, "\n"))
    }
    genofile <- SNPRelate::snpgdsOpen(outfilespec)
    if (verbose >= 3) {
        cat(report("Structure of gds file\n\n"))
        SNPRelate::snpgdsSummary(genofile)
        print(genofile)
    }
    
    # Close the GDS file
    if (verbose >= 2) {
        cat(report("  Closing file", outfilespec, "\n"))
    }
    SNPRelate::snpgdsClose(genofile)
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    invisible(NULL)
    
}
