#' @name gl2hapmap
#' @title Converts a genlight object into HapMap format
#' @family linker
#'
#' @description
#' Convert a \code{genlight} object into HapMap format, producing a file with
#' the standard columns for SNP marker ID, chromosome, position, allele
#' definitions and per-sample genotype calls encoded as nucleotide pairs.
#'
#'  The chromosome information for unmapped SNPs is coded as 0, and the
#'  position of SNPs without chromosome-scale position information is
#'  likewise coded as 0.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file [default 'gl_hapmap'].
#' @param outpath Path where to save the output file [default global working
#' directory or if not specified, tempdir()].
#' @param chrom Field name from the slot loc.metrics where the chromosome of
#' each SNP is stored; if supplied, this field takes precedence over the
#' chromosome slot [default NULL].
#' @param pos Field name from the slot loc.metrics where the SNP position is
#' stored; if supplied, this field takes precedence over the position slot
#' [default NULL].
#' @param strand Orientation of the SNP in the DNA strand. Thus, SNPs could be
#' in the forward (+) or in the reverse (-) orientation relative to the
#' reference genome [default "+"].
#' @param assembly Version of reference sequence assembly (from NCBI)
#' [default ""].
#' @param center Name of genotyping center that produced the genotypes
#' [default NA].
#' @param protLSID Identifier for HapMap protocol [default NA].
#' @param assayLSID Identifier HapMap assay used for genotyping [default NA].
#' @param panelLSID Identifier for panel of individuals genotyped [default NA].
#' @param QCcode Quality control for all entries [default NA].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' \donttest{
#' require("dartR.data")
#' # SNP data
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' gl2hapmap(platypus.gl, outpath = tempdir())
#' }
#'
#' @export
#' @return returns no value (i.e. NULL)

gl2hapmap <- function(x,
                      outfile = "gl_hapmap",
                      outpath = NULL,
                      chrom = NULL,
                      pos = NULL,
                      strand = "+",
                      assembly = "",
                      center = NA,
                      protLSID = NA,
                      assayLSID = NA,
                      panelLSID = NA,
                      QCcode = NA,
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
  # works only with SNP data (F6: the former cat(error()); stop() raised
  # a condition with an empty message)
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  # allele definitions are mandatory for HapMap output (F3)
  if (is.null(x@loc.all)) {
    stop(error(
      "  Fatal Error: the genlight object carries no allele definitions",
      "(the loc.all slot is NULL); allele definitions are required to",
      "write HapMap genotypes.\n"
    ))
  }

  # DO THE JOB
  
  # assigning SNP position information
  # When reading a DArT report, the position of the SNP in the trimmed
  # sequence (presumably always less than 1000) is assigned to the slot
  # position. If the SNP position in the chromosome has been assigned
  # before, from a reference genome for example (presumably always more
  # than 1000), that information is used directly -- unless a loc.metrics
  # field is explicitly nominated via pos, which takes precedence (F1).
  if (!is.null(pos)) {
    metrics <- x$other$loc.metrics
    # field must exist in loc.metrics
    if (!pos %in% names(metrics)) {
      stop(error(sprintf(
        "The field '%s' with SNP position information is not present in loc.metrics.\n",
        pos
      )))
    }
    if (verbose >= 2) {
      cat(report(
        "  Using SNP positions from loc.metrics field '", pos, "'.\n"
      ))
    }
    # pull it out and coerce to integer
    x$position <- as.integer(metrics[[pos]])
  } else if (is.null(x$position) ||
             max(x$position, na.rm = TRUE) < 1000L) {
    # positions absent or within-tag offsets only: zero-fill (F4)
    x$position <- integer(nLoc(x))
    if (verbose >= 2) {
      cat(report(
        "  No chromosome-scale SNP positions available; setting all SNP",
        "positions to '0'.\n"
      ))
    }
  }

  # assigning chromosome information: an explicitly nominated loc.metrics
  # field takes precedence over the chromosome slot (F1)
  if (!is.null(chrom)) {
    metrics <- x$other$loc.metrics
    # require that the chosen field exists
    if (!chrom %in% names(metrics)) {
      stop(error(sprintf(
        "The field '%s' with chromosome information is not present in loc.metrics.\n",
        chrom
      )))
    }
    if (verbose >= 2) {
      cat(report(
        "  Using chromosome data from loc.metrics field '", chrom, "'.\n"
      ))
    }
    # extract and coerce to factor
    x$chromosome <- factor(metrics[[chrom]])
  } else if (is.null(x$chromosome)) {
    # no chromosome field: set all to "0"
    x$chromosome <- factor(rep("0", nLoc(x)))
    if (verbose >= 2) {
      cat(report(
        "  Chromosome slot was NULL; setting all SNP chromosomes to '0'.\n"
      ))
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
  
  x_mat <- as.matrix(x[, ])
  homs1 <-
    paste(substr(x@loc.all, 1, 1), "/", substr(x@loc.all, 1, 1), sep = "")
  hets <- x@loc.all
  homs2 <-
    paste(substr(x@loc.all, 3, 3), "/", substr(x@loc.all, 3, 3), sep = "")
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
        xx[i, ii] <- "N/N"
      }
    }
  }
  xx <- gsub("/", "", xx)
  xx <- as.data.frame(xx)
  xx <- t(xx)
  colnames(xx) <- indNames(x)
  
  geno_tmp <- data.frame(
    rs = locNames(x),
    alleles = x$loc.all,
    chrom = x$chromosome,
    pos = x$position,
    strand = strand,
    assembly = assembly,
    center = center,
    protLSID = protLSID,
    assayLSID = assayLSID,
    panelLSID = panelLSID,
    QCcode = QCcode
  )
  
  colnames(geno_tmp) <- c(
    "rs#"	,
    "alleles",
    "chrom",
    "pos"	,
    "strand",
    "assembly#",
    "center"	,
    "protLSID",
    "assayLSID",
    "panelLSID"	,
    "QCcode"
  )
  
  res_output <- cbind(geno_tmp, xx)
  
  res_output <- res_output[order(res_output$chrom, res_output$pos, decreasing = F), ]
  
  filename1 <- file.path(outpath, paste0(outfile, ".hmp.txt"))
  
  write.table(
    res_output,
    file = filename1,
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = T,
    eol = "\n",
    fileEncoding = "UTF-8"
  )
  
  if (verbose >= 2) {
    cat(report(
      "  The hapmap file is saved as: ",
      filename1, "\n"
    ))
  }

  # FLAG SCRIPT END

  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN
  invisible(NULL)

}
