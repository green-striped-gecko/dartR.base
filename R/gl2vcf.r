#' @name gl2vcf
#' @title Converts a genlight object into vcf format
#' @family linker

#' @description
#' This function exports a genlight object into VCF format and save it into a
#' file.
#'
#' @details
#' This function requires the binary file of PLINK 1.9 to be downloaded and
#' its path provided (plink.bin.path).
#' The binary file can be downloaded from:
#' \url{https://www.cog-genomics.org/plink/}
#'
#' The slots \code{@position} and \code{@chromosome} hold genome
#' coordinates. Explicitly supplied snp.pos/snp.chr loc.metrics fields take
#' precedence over the slots; when neither a valid slot nor a field is
#' available, CHROM and POS are written as 0 (unmapped).
#'
#' Family ID is taken from  x$pop
#'
#' Within-family ID (cannot be '0') is taken from indNames(x)
#'
#' Variant identifier is taken from locNames(x)
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param plink.bin.path Path of PLINK binary file [default getwd()].
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
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @return  returns no value (i.e. NULL)
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
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
  # SNP only: the ped/map recoding requires loc.all and 0/1/2 dosages,
  # neither of which SilicoDArT (presence/absence) data carries
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  # PLINK 1.9 must be present before any work is done
  if (!any(file.exists(file.path(plink.bin.path, c("plink", "plink.exe"))))) {
    stop(error(
      "Cannot find the PLINK executable ('plink' or 'plink.exe') in '",
      plink.bin.path, "'. This function needs PLINK 1.9 to work. Please ",
      "download it from https://www.cog-genomics.org/plink/ and supply ",
      "its directory via plink.bin.path.\n"
    ))
  }

  # DO THE JOB

  # assigning SNP position and chromosome information
  # The @position/@chromosome slots are reserved for genome coordinates
  # and are NULL until assigned (the position of the SNP within the
  # sequence tag lives in @other$loc.metrics$SnpPosition). An explicitly
  # supplied snp.pos/snp.chr field takes precedence over the slots;
  # otherwise a valid slot is used as found; otherwise coordinates are
  # written as 0 (unmapped).
  metrics <- x$other$loc.metrics

  if (!is.null(snp.pos)) {
    # field must exist in loc.metrics
    if (!snp.pos %in% names(metrics)) {
      stop(error(sprintf(
        "The field '%s' with SNP position information is not present in loc.metrics.\n",
        snp.pos
      )))
    }
    # coerce via character so a factor-typed field yields its labels,
    # not its level codes; refuse values that are not integer positions
    pos.tmp <- suppressWarnings(as.integer(as.character(metrics[[snp.pos]])))
    if (any(is.na(pos.tmp))) {
      stop(error(sprintf(
        "The field '%s' contains values that cannot be interpreted as integer SNP positions.\n",
        snp.pos
      )))
    }
    if (verbose >= 1) {
      cat(report(
        "  Using SNP positions from loc.metrics field '", snp.pos, "'\n"
      ))
    }
    x$position <- pos.tmp
  } else if (!is.null(x$position) && length(x$position) == nLoc(x)) {
    # a valid slot is used as found
    if (verbose >= 1) {
      cat(report("  Using SNP positions from the @position slot\n"))
    }
  } else {
    # no usable source: zero out (unmapped), warning if a malformed
    # slot is being discarded
    if (!is.null(x$position) && verbose >= 1) {
      cat(warn(
        "  Warning: @position slot has length", length(x$position),
        "but the object has", nLoc(x), "loci; positions set to 0\n"
      ))
    }
    x$position <- integer(nLoc(x))
  }

  if (!is.null(snp.chr)) {
    # require that the chosen field exists
    if (!snp.chr %in% names(metrics)) {
      stop(error(sprintf(
        "The field '%s' with chromosome information is not present in loc.metrics.\n",
        snp.chr
      )))
    }
    if (verbose >= 1) {
      cat(report(
        "  Using chromosome data from loc.metrics field '", snp.chr, "'\n"
      ))
    }
    x$chromosome <- factor(metrics[[snp.chr]])
  } else if (!is.null(x$chromosome) && length(x$chromosome) == nLoc(x)) {
    # a valid slot is used as found
    if (verbose >= 1) {
      cat(report("  Using chromosomes from the @chromosome slot\n"))
    }
  } else {
    # no usable source: set all to "0" (unmapped), warning if a
    # malformed slot is being discarded
    if (!is.null(x$chromosome) && verbose >= 1) {
      cat(warn(
        "  Warning: @chromosome slot has length", length(x$chromosome),
        "but the object has", nLoc(x), "loci; chromosomes set to '0'\n"
      ))
    }
    x$chromosome <- factor(rep("0", nLoc(x)))
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
  
  # the ped/map intermediates are not part of the requested output:
  # write them to the per-session tempdir() rather than outpath
  gl2plink(
    x = x,
    outfile = "gl_plink_temp",
    outpath = tempdir(),
    chr.format = chr.format,
    pos.cM = pos.cM,
    ID.dad = ID.dad,
    ID.mum = ID.mum,
    sex.code = sex.code,
    phen.value = phen.value,
    verbose = verbose
  )

  prefix.in_temp <- file.path(tempdir(), "gl_plink_temp")
  prefix.out_temp <- file.path(outpath, outfile)

  # allele list for PLINK: column 2 = alternate allele (forced A1),
  # column 3 = reference allele (forced A2), so REF/ALT always come from
  # loc.all, including loci fixed for the alternate allele in the sample
  allele_tmp <- gsub("/", " ", x$loc.all)
  allele_tmp <- strsplit(allele_tmp, split = " ")
  allele_tmp <- Reduce(rbind, allele_tmp)
  allele_tmp <- cbind(locNames(x), allele_tmp[, 2], allele_tmp[, 1])
  write.table(
    allele_tmp,
    file = file.path(tempdir(), "mylist.txt"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  # PLINK 1.9 rejects --a1-allele (--reference-allele) and --a2-allele in
  # the same run, so the alleles are pinned in two passes: A1 (= ALT)
  # while converting to bed, then A2 (= REF) while recoding to VCF. This
  # keeps REF/ALT faithful to loc.all even for loci fixed for one allele
  # in the exported sample.
  prefix.bed_temp <- file.path(tempdir(), "gl_plink_temp_bed")
  make_plink <-
    function(plink.bin.path,
             prefix.in = prefix.in_temp,
             prefix.out = prefix.out_temp,
             autosome.only = FALSE,
             extra.options = "") {
      system_verbose(
        paste(
          plink.bin.path,
          "--file",
          prefix.in,
          "--make-bed",
          if (autosome.only)
            "--autosome"
          else
            "",
          "--allow-no-sex",
          paste("--a1-allele", file.path(tempdir(), 'mylist.txt'), "2", "1"),
          "--out",
          prefix.bed_temp,
          extra.options
        )
      )
      system_verbose(
        paste(
          plink.bin.path,
          "--bfile",
          prefix.bed_temp,
          "--recode",
          "vcf",
          "--allow-no-sex",
          paste("--a2-allele", file.path(tempdir(), 'mylist.txt'), "3", "1"),
          "--out",
          prefix.out,
          extra.options
        )
      )
    }

  # PLINK's console log is progress information: keep it captured for
  # error diagnosis but print it only at verbose >= 2 (PLINK's stderr
  # stream is likewise suppressed below that level)
  system_verbose <- function(...) {
    plink.log <- system(..., intern = TRUE, ignore.stderr = (verbose < 2))
    if (verbose >= 2) {
      message(
        paste0(
          "\n\n----------Output of function start:\n\n",
          paste(plink.log, collapse = "\n"),
          "\n\n----------Output of function finished...\n\n"
        )
      )
    }
    invisible(plink.log)
  }

  make_plink(plink.bin.path = paste0(plink.bin.path, "/plink"),
             extra.options = "--aec")

  # tidy up: on success, remove the ped/map/bed intermediates and the
  # PLINK by-products (.log, .nosex); on failure they are left for
  # diagnosis
  if (file.exists(paste0(prefix.out_temp, ".vcf"))) {
    unlink(paste0(prefix.in_temp, c(".map", ".ped")))
    unlink(paste0(prefix.bed_temp,
                  c(".bed", ".bim", ".fam", ".log", ".nosex")))
    unlink(paste0(prefix.out_temp, c(".log", ".nosex")))
  }
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  invisible(NULL)
}
