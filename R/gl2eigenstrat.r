#' @name gl2eigenstrat
#' @title Converts a genlight object into eigenstrat format
#' @family linker

#' @description
#' The output of this function are three files:
#' \itemize{
#' \item genotype file: contains genotype data for each individual at each SNP
#' with an extension 'eigenstratgeno.'
#' \item snp file: contains information about each SNP with an extension 'snp.'
#' \item indiv file: contains information about each individual with an
#' extension 'ind.'
#' }
#' 
#' @details
#' Each value in the genotype file is the number of copies of the reference
#' allele -- the first allele of \code{loc.all}, written in the fifth column
#' of the snp file: 2 homozygous reference, 1 heterozygous, 0 homozygous
#' alternate, 9 missing.
#'
#' Eigenstrat only accepts chromosomes coded as numeric values, as follows:
#' X chromosome is encoded as 23, Y is encoded as 24, mtDNA is encoded as
#' 90, and XY is encoded as 91. Chromosome labels 'X', 'Y', 'MT'/'mtDNA'
#' and 'XY' in the nominated field are mapped to this encoding; SNPs with
#' illegal chromosome values (any other non-numeric label, or a value less
#' than 1, such as 0) are removed, with a warning stating how many. If no
#' locus has an encodable chromosome value, the function stops.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file [default 'gl_eigenstrat'].
#' @param outpath Path where to save the output file [default global working
#' directory or if not specified, tempdir()].
#' @param snp.pos Field name from the slot loc.metrics where the SNP position
#' is stored; when left at the numeric sentinel 1, a constant position of 1
#' is written for every SNP [default 1].
#' @param snp.chr Field name from the slot loc.metrics where the chromosome
#' of each SNP is stored; when left at the numeric sentinel 1, a constant
#' chromosome of 1 is written for every SNP [default 1].
#' @param pos.cM A vector, with as many elements as there are loci, containing
#' the SNP position in morgans or centimorgans; the default writes 0 (unknown)
#' for every SNP [default 0].
#' @param sex.code A vector, with as many elements as there are individuals,
#' containing the sex code ('male', 'female', 'unknown') [default 'unknown'].
#' @param phen.value A vector, with as many elements as there are individuals,
#' containing the phenotype value ('Case', 'Control') [default 'Case'].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @return returns no value (i.e. NULL)
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @references
#' \itemize{
#' \item Patterson, N., Price, A. L., & Reich, D. (2006). Population structure
#' and eigenanalysis. PLoS genetics, 2(12), e190.
#' \item Price, A. L., Patterson, N. J., Plenge, R. M., Weinblatt, M. E.,
#' Shadick, N. A., & Reich, D. (2006). Principal components analysis corrects
#' for stratification in genome-wide association studies. Nature genetics,
#' 38(8), 904-909.
#' }
#'
#' @examples
#' \donttest{
#' require("dartR.data")
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' # The platypus chromosome labels (e.g. 'NC_041731.1_chromosome_4') cannot
#' # be encoded for Eigenstrat, so the chromosome field is left at its
#' # constant default
#' gl2eigenstrat(platypus.gl,snp.pos='ChromPos_Platypus_Chrom_NCBIv1',
#' outpath=tempdir())
#' }
#'
#' @export

gl2eigenstrat <- function(x,
                          outfile = "gl_eigenstrat",
                          outpath = NULL,
                          snp.pos = 1,
                          snp.chr = 1,
                          pos.cM = 0,
                          sex.code = "unknown",
                          phen.value = "Case",
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

    # FUNCTION SPECIFIC ERROR CHECKING

    # cbind would silently recycle short sex.code/phen.value vectors down
    # the .ind file [approved F6]
    if (!length(sex.code) %in% c(1, nInd(x))) {
        stop(error(
            "Fatal Error: sex.code must be a single value or a vector with",
            "one element per individual\n"
        ))
    }
    if (!length(phen.value) %in% c(1, nInd(x))) {
        stop(error(
            "Fatal Error: phen.value must be a single value or a vector with",
            "one element per individual\n"
        ))
    }

    # DO THE JOB

    # Resolve the chromosome of each SNP first, so that loci with illegal
    # chromosome values can be removed before any file is written
    snp_temp <- x$other$loc.metrics

    if (snp.chr == 1) {
        chrom <- 1
    } else {
        # Coerce via as.character so a factor field yields its labels, not
        # its level codes [approved F2], then map the documented non-numeric
        # encodings
        chrom_raw <- as.character(unname(unlist(snp_temp[snp.chr])))
        chrom_raw[toupper(chrom_raw) == "X"] <- "23"
        chrom_raw[toupper(chrom_raw) == "Y"] <- "24"
        chrom_raw[toupper(chrom_raw) %in% c("MT", "MTDNA")] <- "90"
        chrom_raw[toupper(chrom_raw) == "XY"] <- "91"
        chrom <- suppressWarnings(as.numeric(chrom_raw))

        # Remove SNPs with illegal chromosome values (non-numeric after the
        # mapping, or less than 1) as documented [approved F3]; if no locus
        # can be encoded, stop rather than write an empty file set
        # [approved F2]
        illegal <- is.na(chrom) | chrom < 1
        if (all(illegal)) {
            stop(error(
                "Fatal Error: no locus has a chromosome value that can be",
                "encoded for Eigenstrat; the nominated field must hold",
                "numeric labels or 'X', 'Y', 'MT'/'mtDNA', 'XY'\n"
            ))
        }
        if (any(illegal)) {
            if (verbose >= 1) {
                cat(warn(
                    "  Warning:",
                    sum(illegal),
                    "loci with illegal chromosome values removed from the",
                    "output\n"
                ))
            }
            keep <- which(!illegal)
            old.metrics <- x@other$loc.metrics
            x <- x[, keep]
            x@other$loc.metrics <- old.metrics[keep, , drop = FALSE]
            snp_temp <- x@other$loc.metrics
            chrom <- chrom[keep]
            if (length(pos.cM) > 1) {
                pos.cM <- pos.cM[keep]
            }
        }
    }

    if (snp.pos == 1) {
        snp.pos <- 1
    } else {
        # as.character first, for the same factor-coercion reason as the
        # chromosome field [approved F2]
        snp.pos <- as.numeric(as.character(unname(unlist(snp_temp[snp.pos]))))
    }

    # geno file
    geno_temp <- t(as.matrix(x))
    # EIGENSTRAT geno values count copies of the REFERENCE allele (the first
    # loc.all allele, written in .snp column 5), while the dartR score counts
    # the alternate allele -- invert the score [approved F1]
    geno_temp <- 2 - geno_temp
    geno_temp[is.na(geno_temp)] <- 9
    geno_file <-
        as.matrix(unname(unlist(
            apply(geno_temp, 1, paste0, collapse = "")
        )))

    write.table(
        geno_file,
        file = paste0(outfilespec, ".eigenstratgeno"),
        quote = F,
        row.names = F,
        col.names = F
    )
    # snp file
    snp_name <- locNames(x)
    
    ref_allele <- substring(x@loc.all, 1, 1)
    var_allele <- substring(x@loc.all, 3, 3)
    
    snp_file <-
        cbind(snp_name, chrom, pos.cM, snp.pos, ref_allele, var_allele)
    
    write.table(
        snp_file,
        file = paste0(outfilespec, ".snp"),
        quote = F,
        row.names = F,
        col.names = F
    )
    
    # indiv file
    sample_id <- indNames(x)
    
    # 2nd column is gender (M or F).  If unknown, ok to set to U for Unknown.
    sex.code[sex.code == "female"] <- "F"
    sex.code[sex.code == "male"] <- "M"
    sex.code[sex.code == "unknown"] <- "U"
    
    indiv_file <- cbind(sample_id, sex.code, phen.value)
    
    write.table(
        indiv_file,
        file = paste0(outfilespec, ".ind"),
        quote = F,
        row.names = F,
        col.names = F
    )
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    invisible(NULL)
}
