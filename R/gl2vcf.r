#' @name gl2vcf
#' @title Converts a genlight object into vcf format
#' @family linker

#' @description
#' This function exports a genlight object into VCF format and save it into a
#' file.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param plink.bin.path Path of PLINK binary file [default getwd())].
#' @param outfile File name of the output file [default 'gl_vcf'].
#' @param outpath Path where to save the output file [default global working
#' directory or if not specified, tempdir()].
#' @param snp.pos Field name from the slot loc.metrics where the SNP position is
#' stored [default NULL].
#' @param snp.chr Field name from the slot loc.metrics where the chromosome of
#' each is stored [default NULL].
#' @param chr.format Whether chromosome information is stored as 'numeric' or as
#' 'character', see details [default 'character'].
#' @param pos.cM A vector, with as many elements as there are loci, containing
#' the SNP position in morgans or centimorgans [default '0'].
#' @param ID.dad A vector, with as many elements as there are individuals,
#' containing the ID of the father, '0' if father isn't in dataset [default '0'].
#' @param ID.mum A vector, with as many elements as there are individuals,
#' containing the ID of the mother, '0' if mother isn't in dataset [default '0'].
#' @param sex.code A vector, with as many elements as there are individuals,
#' containing the sex code ('male', 'female', 'unknown') [default  'unknown'].
#' @param phen.value A vector, with as many elements as there are individuals,
#' containing the phenotype value. '1' = control, '2' = case, '0' = unknown
#' [default '0'].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#'
#' @details
#' This function requires to download the binary file of PLINK 1.9 and provide
#' its path (plink.bin.path).
#' The binary file can be downloaded from:
#' \url{https://www.cog-genomics.org/plink/}
#'
#' The chromosome information for unmapped SNPS is coded as 0.
#'
#' Family ID is taken from  x$pop
#'
#' Within-family ID (cannot be '0') is taken from indNames(x)
#'
#' Variant identifier is taken from locNames(x)
#'
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' \dontrun{
#' # this example needs plink installed to work
#' require("dartR.data")
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' gl2vcf(platypus.gl,snp.pos='ChromPos_Platypus_Chrom_NCBIv1',
#'  snp.chr = 'Chrom_Platypus_Chrom_NCBIv1')
#' }
#'
#' @references
#' Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M.
#'  A., ... & 1000 Genomes Project Analysis Group. (2011). The variant call
#'  format and VCFtools. Bioinformatics, 27(15), 2156-2158.
#'
#' @export
#' @return  returns no value (i.e. NULL)

gl2vcf <- function(x,
                   plink.bin.path = getwd(),
                   outfile = "gl_vcf",
                   outpath = NULL,
                   snp.pos = NULL,
                   snp.chr = NULL,
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
  outpath <- gl.check.wd(outpath, verbose = 0)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # DO THE JOB
  
  # assigning SNP position information
  # When reading a DArT report, the position of the SNP in the trimmedsequence
  # (presumably is always less than 1000) is assigned to the slot position.
  # If the SNP position in the chromosome has been assigned before from a
  # reference genome for example  (presumably always more than 1000), use that
  # information directly.
  
  # only reset positions if the current max is < 1000
  if (max(x$position, na.rm = TRUE) < 1000L) {
    # no SNP‐position field supplied: zero out
    if (is.null(snp.pos)) {
      x$position <- integer(nLoc(x))
    } else {
      metrics <- x$other$loc.metrics
      # field must exist in loc.metrics
      if (!snp.pos %in% names(metrics)) {
        stop(error(sprintf(
          "The field '%s' with SNP position information is not present in loc.metrics.\n",
          snp.pos
        )))
      }
      # verbose message
      if (verbose >= 2) {
        message(report(
          "Using SNP positions from loc.metrics field '", snp.pos, "'.\n"
        ))
      }
      # pull it out and coerce to integer
      x$position <- as.integer(metrics[[snp.pos]])
    }
  }
  
  # assign chromosome information if missing
  if (is.null(x$chromosome)) {
    metrics <- x$other$loc.metrics
    if (is.null(snp.chr)) {
      # no chromosome field: set all to "0"
      x$chromosome <- factor(rep("0", nLoc(x)))
      if (verbose >= 2) {
        message(report(
          "Chromosome slot was NULL; setting all SNP chromosomes to '0'.\n"
        ))
      }
    } else {
      # require that the chosen field exists
      if (!snp.chr %in% names(metrics)) {
        stop(error(sprintf(
          "The field '%s' with chromosome information is not present in loc.metrics.\n",
          snp.chr
        )))
      }
      if (verbose >= 2) {
        message(report(
          "Using chromosome data from loc.metrics field '", snp.chr, "'.\n"
        ))
      }
      # extract and coerce to factor
      x$chromosome <- factor(metrics[[snp.chr]])
    }
  }
  
  # Chromosome "0" is assigned to unmmapped SNPs
  # ensure "0" is a valid level
  if (!"0" %in% levels(x$chromosome)) {
    levels(x$chromosome) <- c(levels(x$chromosome), "0")
  }
  #replace blanks
  x$chromosome[x$chromosome == ""] <- "0"
  #drop any now‐unused levels
  x$chromosome <- droplevels(x$chromosome)
  
  gl2plink(
    x = x,
    outfile = "gl_plink_temp",
    outpath = outpath,
    chr.format = chr.format,
    pos.cM = pos.cM,
    ID.dad = ID.dad,
    ID.mum = ID.mum,
    sex.code = sex.code,
    phen.value = phen.value,
    verbose = NULL
  )
  
  prefix.in_temp <- paste0(outpath, "/gl_plink_temp")
  prefix.out_temp <- file.path(outpath, outfile)
  
  allele_tmp <- gsub("/", " ", x$loc.all)
  allele_tmp <- strsplit(allele_tmp, split = " ")
  allele_tmp <- Reduce(rbind, allele_tmp)[, 2]
  allele_tmp <- cbind(locNames(x), allele_tmp)
  write.table(
    allele_tmp,
    file = file.path(tempdir(), "mylist.txt"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  make_plink <-
    function(plink.bin.path,
             prefix.in = prefix.in_temp,
             prefix.out = prefix.out_temp,
             autosome.only = FALSE,
             extra.options = "") {
      bedfile.out <- paste0(prefix.out, ".bed")
      system_verbose(
        paste(
          plink.bin.path,
          "--file",
          prefix.in,
          "--recode",
          "vcf",
          if (autosome.only)
            "--autosome"
          else
            "",
          "--allow-no-sex",
          paste("--reference-allele", file.path(tempdir(), 'mylist.txt')),
          # "--keep-allele-order",
          # "--real-ref-alleles",
          # paste("--a1-allele", file.path(outpath,'alleles.csv'),"1"),
          # paste("--a2-allele", file.path(outpath,'alleles.csv'),"2"),
          "--out",
          prefix.out,
          extra.options
        )
      )
    }
  
  system_verbose <- function(...) {
    report <- system(..., intern = T)
    message(
      paste0(
        "\n\n----------Output of function start:\n\n",
        paste(report, collapse = "\n"),
        "\n\n----------Output of function finished...\n\n"
      )
    )
  }
  
  make_plink(plink.bin.path = paste0(plink.bin.path, "/plink"),
             extra.options = "--aec")
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  invisible(NULL)
}
