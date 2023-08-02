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
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file [default 'gl_eigenstrat'].
#' @param outpath Path where to save the output file [default global working 
#' directory or if not specified, tempdir()].
#' @param snp.pos Field name from the slot loc.metrics where the SNP position is
#' stored [default 1].
#' @param snp.chr Field name from the slot loc.metrics where the chromosome of
#' each is stored [default 1].
#' @param pos.cM A vector, with as many elements as there are loci, containing
#' the SNP position in morgans or centimorgans [default 1].
#' @param sex.code A vector, with as many elements as there are individuals,
#' containing the sex code ('male', 'female', 'unknown') [default 'unknown'].
#' @param phen.value A vector, with as many elements as there are individuals,
#' containing the phenotype value ('Case', 'Control') [default 'Case'].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' 
#' @details
#' Eigenstrat only accepts chromosomes coded as numeric values, as follows:
#' X chromosome is encoded as 23, Y is encoded as 24, mtDNA is encoded as
#' 90, and XY is encoded as 91. SNPs with illegal chromosome values, such
#' as 0, will be removed.
#' 
#' @author Custodian: Luis Mijangos 
#' (Post to \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#' \donttest{
#' require("dartR.data")
#' gl2eigenstrat(platypus.gl,snp.pos='ChromPos_Platypus_Chrom_NCBIv1',
#' snp.chr = 'Chrom_Platypus_Chrom_NCBIv1', outpath=tempdir())
#' }
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
#' @return  returns no value (i.e. NULL)
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
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # DO THE JOB
    
    # geno file
    geno_temp <- t(as.matrix(x))
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
    
    snp_temp <- x$other$loc.metrics
    
    if (snp.chr == 1) {
        chrom <- 1
    } else {
        chrom <- as.numeric(unname(unlist(snp_temp[snp.chr])))
    }
    
    if (snp.pos == 1) {
        snp.pos <- 1
    } else {
        snp.pos <- as.numeric(unname(unlist(snp_temp[snp.pos])))
    }
    
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
