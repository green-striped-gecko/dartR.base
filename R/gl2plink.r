#' @name gl2plink
#' @title Converts a genlight object into PLINK format
#' @family linkers

#' @description
#' This function exports a genlight object into PLINK format and save it into a
#' file.
#' By default it produces two files: a .ped file (which embeds the six sample
#' columns of the .fam format) and a .map file. With bed.files = TRUE and an
#' external PLINK 1.9 binary, the binary .bed, .bim and .fam files are also
#' produced.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param plink.bin.path Path of PLINK binary file [default getwd()].
#' @param bed.files Whether create PLINK files .bed, .bim and .fam
#' [default FALSE].
#' @param outfile File name of the output file [default 'gl_plink'].
#' @param outpath Path where to save the output file [default global working 
#' directory or if not specified, tempdir()].
#' @param chr.format Whether chromosome information is stored as 'numeric' or as
#' 'character', see details [default 'character'].
#' @param pos.cM A vector, with as many elements as there are loci, containing
#' the SNP position in morgans or centimorgans [default '0'].
#' @param ID.dad A vector, with as many elements as there are individuals,
#' containing the ID of the father, '0' if father isn't in dataset [default '0'].
#' @param ID.mum A vector, with as many elements as there are individuals,
#' containing the ID of the mother, '0' if mother isn't in dataset [default '0'].
#' @param sex.code A vector, with as many elements as there are individuals,
#' containing the sex code ('male', 'female', 'unknown'). Sex information needs 
#' just to start with an "F" or "f" for females, with an "M" or "m" for males 
#' and with a "U", "u" or being empty if the sex is unknown [default 'unknown'].
#' @param phen.value A vector, with as many elements as there are individuals,
#' containing the phenotype value. '1' = control, '2' = case, '0' = unknown
#' [default '0'].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @details
#' To create PLINK files .bed, .bim and .fam (bed.files = TRUE), it is necessary
#' to download the binary file of PLINK 1.9 and provide its path (plink.bin.path).
#' The binary file can be downloaded from:
#' \url{https://www.cog-genomics.org/plink/}
#' 
#' After downloading, unzip the file, access the unzipped folder and move the 
#' binary file ("plink") to your working directory.
#' 
#' If you are using a Mac, you might need to open the binary first to grant 
#' access to the binary. 
#' 
#' The chromosome of each SNP can be a character or numeric. The chromosome
#' information for unmapped SNPS is coded as 0.
#'
#' If the accessors x$chromosome and x$position are unset (NULL), the loci are
#' treated as unmapped: chromosome and base-pair position are written as 0,
#' PLINK's own convention, and a warning is given. Note also that if
#' x$position holds within-tag offsets (as when read directly from a DArT
#' report, values typically < 70) rather than genomic coordinates, those
#' offsets are written to the .map as base-pair positions; assign genomic
#' coordinates before calling if downstream position-aware PLINK operations
#' are intended.
#'
#' Family ID is taken from  x$pop.
#' 
#' Within-family ID (cannot be '0') is taken from indNames(x).
#' 
#' Variant identifier is taken from locNames(x).
#' 
#' SNP position is taken from the accessor x$position.
#' 
#' Chromosome name is taken from the accessor x$chromosome
#' 
#' Note that if names of populations or individuals contain spaces, they are 
#' replaced by an underscore "_".
#' 
#' If you like to use chromosome information when converting to plink format and
#'  your chromosome names are not from human, you need to change the chromosome 
#'  names as 'contig1', 'contig2', etc. as described in the section "Nonstandard
#'   chromosome IDs" in the following link:
#'   https://www.cog-genomics.org/plink/1.9/input
#'
#' @return returns no value (i.e. NULL)
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#' require("dartR.data")
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' test <- platypus.gl
#' # assigning SNP position
#' test$position <- test$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
#' # assigning a dummy name for chromosomes
#' test$chromosome <- as.factor("1")
#' gl2plink(test, outpath=tempdir())
#' }
#' 
#' @references
#' Purcell, Shaun, et al. 'PLINK: a tool set for whole-genome association and
#' population-based linkage analyses.' The American journal of human genetics
#' 81.3 (2007): 559-575.
#'
#' @export

gl2plink <- function(x,
                     plink.bin.path = getwd(),
                     bed.files = FALSE,
                     outfile = "gl_plink",
                     outpath = NULL,
                     chr.format = "character",
                     pos.cM = "0",
                     ID.dad = "0",
                     ID.mum = "0",
                     sex.code = "unknown",
                     phen.value = "0",
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
    # SNP only: the recoding requires loc.all and 0/1/2 dosage (F3)
    datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

    # DO THE JOB

    # Unset slots are treated as unmapped loci and written with PLINK's own
    # convention, chromosome 0 and position 0, instead of fabricated
    # coordinates (F1); results-affecting, so warned at verbose >= 1 (F5)
    if(is.null(x$chromosome)){
      x$chromosome <- as.factor(rep("0",nLoc(x)))
      if (verbose >= 1) {
        cat(warn("  Warning: chromosome slot is empty. Loci are coded as unmapped (chromosome 0).\n"))
      }
    }

    if(is.null(x$position)){
      x$position <- rep(0L, nLoc(x))
      if (verbose >= 1) {
        cat(warn("  Warning: position slot is empty. Loci are coded as unmapped (position 0).\n"))
      }
    }
    
    snp_temp <- as.data.frame(cbind(as.character(x$chromosome),
                                    as.numeric(x$position)))
    colnames(snp_temp) <- c("chrom","snp_pos")
    

        if (chr.format == "numeric") {
            snp_temp$chrom <- as.numeric(snp_temp$chrom)
        }
        if (chr.format == "character") {
            snp_temp$chrom <- as.character(snp_temp$chrom)
        }
    
        snp_temp$snp_pos <- as.numeric(snp_temp$snp_pos)
    
    # Convert any NA values to 0
    snp_temp[is.na(snp_temp$snp_pos), "snp_pos"] <- 0
    # Convert any NA values to 0
    snp_temp[snp_temp$snp_pos == 0, "chrom"] <- 0
    
    # Chromosome code
    snp_chr <- snp_temp$chrom
    # Variant identifier
    var_id <- locNames(x)
    
    # Base-pair coordinate
    pos_bp <- snp_temp$snp_pos
    
    gl_map <- cbind(snp_chr, var_id, pos.cM, pos_bp)
    
    write.table(
        gl_map,
        file = paste0(outfilespec, ".map"),
        quote = F,
        row.names = F,
        col.names = F
    )
    
    ########## .fam (PLINK sample information file)
    
    sample.id_temp <- indNames(x)
    sample.id_temp <-
        gsub(" ", replacement = "_", sample.id_temp)
    
    # Family ID ('FID')
    FID <- as.character(x$pop)
    FID <- gsub(" ","_",FID)
    # Within-family ID ('IID'; cannot be '0')
    IID <- sample.id_temp

    ID.dad <- gsub(" ","_",ID.dad)
    ID.mum <- gsub(" ","_",ID.mum)

    # Sex code ('1' = male, '2' = female, '0' = unknown)
    # Recoded unconditionally; scalars (including the default 'unknown')
    # previously reached the .ped verbatim (F2)
    sex.code <- as.character(sex.code)
    sex.code[startsWith(sex.code, "f") |
                 startsWith(sex.code, "F")] <- "2"
    sex.code[startsWith(sex.code, "m") |
                 startsWith(sex.code, "M")] <- "1"
    sex.code[startsWith(sex.code, "u") |
                 startsWith(sex.code, "U")] <- "0"
    sex.code[nchar(sex.code) == 0] <- "0"
    
    gl_fam <-
        cbind(FID, IID, ID.dad, ID.mum, sex.code, phen.value)
    
    x_mat <- as.matrix(x[, ])
    homs1 <- paste(substr(x@loc.all, 1, 1), "/", substr(x@loc.all, 1, 1), sep = "")
    hets <- x@loc.all
    homs2 <- paste(substr(x@loc.all, 3, 3), "/", substr(x@loc.all, 3, 3), sep = "")
    xx <- matrix(NA, ncol = ncol(x_mat), nrow = nrow(x_mat))
    for (i in 1:nrow(x_mat)) {
        for (ii in 1:ncol(x_mat)) {
            inp <- x_mat[i, ii]
            if (!is.na(inp)) {
                if (inp == 0)
                    xx[i, ii] <- homs1[ii]
                else if (inp == 1)
                    xx[i, ii] <- hets[ii]
                else if (inp == 2)
                    xx[i, ii] <- homs2[ii]
            } else{
                xx[i, ii] <-"0/0"
            }
        }
    }
    xx <- gsub("/", " ", xx)
    xx <- apply(xx, 1, paste, collapse = " ")
    xx <- cbind(gl_fam, xx)
    
    write.table(
        xx,
        file = paste0(outfilespec, ".ped"),
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )
    
    if (bed.files) {
        prefix.in_temp <- outfilespec
        prefix.out_temp <- outfilespec
        
        # Reference (A2) allele list: locus name, first allele of loc.all.
        # Named from outfile and passed to PLINK via --a2-allele so the
        # .bed/.bim ref/alt orientation matches the genlight, rather than
        # being reassigned by minor-allele frequency (F4)
        allele_tmp <- gsub("/"," ", x$loc.all)
        allele_tmp <- strsplit(allele_tmp,split = " ")
        allele_tmp <- Reduce(rbind,allele_tmp)[,1]
        allele_tmp <- cbind(locNames(x), allele_tmp)
        allele.file <- file.path(outpath, paste0(outfile, "_a2_alleles.txt"))
        write.table(allele_tmp,
                    file = allele.file,
                    row.names = FALSE,
                    col.names = FALSE,
                    quote = FALSE
                    )

        make_plink <-
            function(plink.path,
                     prefix.in = prefix.in_temp,
                     prefix.out = prefix.out_temp,
                     extra.options = "") {
                bedfile.out <- paste0(prefix.out, ".bed")
                system_verbose(
                    paste(
                        plink.path,
                        "--file",
                        prefix.in,
                        "--allow-no-sex",
                        "--allow-extra-chr",
                        # A2 allele in column 2, variant ID in column 1 (F4)
                        "--a2-allele",
                        allele.file,
                        "2",
                        "1",
                        "--out",
                        prefix.out,
                        extra.options
                    )
                )
                bedfile.out
            }
        
        
        system_verbose <-function(...) {
            report <-system(..., intern = T)
            message(
                paste0(
                    "\n\n----------Output of function start:\n\n",
                    paste(report, collapse = "\n"),
                    "\n\n----------Output of function finished...\n\n"
                )
            )
        }
        
        make_plink(plink.path = paste0(plink.bin.path, "/plink"))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    invisible(NULL)
}
