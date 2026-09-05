#' @name gl2gapit
#' @title Converts a genlight object into a format suitable for input to GAPIT
#' @family linker

#' @description
#' Converts a genlight object containing SNP data into a hapmap-style
#' data.frame as expected by GAPIT (Genome Association and Prediction
#' Integrated Tool) and writes it to file as tab-delimited text.

#' @details
#' The returned data.frame follows the hapmap layout that GAPIT reads with
#' header = FALSE: the first row repeats the column names, the first eleven
#' columns are the hapmap descriptors (rs, alleles, chrom, pos, strand,
#' assembly, center, protLSID, assayLSID, panel, QCcode) and each remaining
#' column carries one individual. Genotypes double the allele letters (for
#' example CC, CT) and missing data is coded 00.
#'
#' Chromosome names are recoded to integer codes assigned in alphabetical
#' order of the distinct names; the name-to-code mapping is reported at
#' verbose >= 2. If the chromosome or position slot is empty, dummy values
#' are used (chromosome 1 for all loci; positions 1 to the number of loci)
#' with a warning at verbose >= 2.
#'
#' The data.frame is written to file.path(outpath, outfile) as tab-delimited
#' text without quotes, ready for GAPIT to read.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file [default 'gl_gapit'].
#' @param outpath Path where to save the output file [default global working
#' directory or if not specified, tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @return The hapmap-style data.frame, returned invisibly.
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#' t1 <- platypus.gl
#' # assigning chromosome
#' t1$chromosome <- t1$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1
#' # assigning SNP position
#' t1$position <- t1$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
#' res <- gl2gapit(t1, outpath = tempdir())
#' }
#'
#' @export

gl2gapit <- function(x,
                     outfile = "gl_gapit",
                     outpath = NULL,
                     verbose = NULL){

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

  x_mat <- as.matrix(x[, ])
  homs1 <- paste(substr(x@loc.all, 1, 1), "/", substr(x@loc.all, 1, 1), sep = "")
  hets <- x@loc.all
  homs2 <- paste(substr(x@loc.all, 3, 3), "/", substr(x@loc.all, 3, 3), sep = "")
  # vectorised genotype-to-letter mapping: index the per-locus codes by
  # column; missing data is coded "0/0"
  xx <- matrix("0/0", nrow = nrow(x_mat), ncol = ncol(x_mat))
  cix <- col(x_mat)
  sel <- !is.na(x_mat) & x_mat == 0
  xx[sel] <- homs1[cix[sel]]
  sel <- !is.na(x_mat) & x_mat == 1
  xx[sel] <- hets[cix[sel]]
  sel <- !is.na(x_mat) & x_mat == 2
  xx[sel] <- homs2[cix[sel]]
  xx <- gsub("/", "", xx)
  xx <- as.data.frame(xx)
  xx <- t(xx)
  colnames(xx) <- indNames(x)

  # If the chromosome/position slots are empty, fill them with dummies as
  # gl2plink does - the hapmap data.frame below needs one value per locus,
  # and empty slots error with "arguments imply differing number of rows".
  if (is.null(x$chromosome) || length(x$chromosome) == 0) {
    x$chromosome <- as.factor(rep("1", nLoc(x)))
    if (verbose >= 2) {
      cat(warn("  Chromosome slot is empty. Using 1 as dummy name.\n"))
    }
  }
  if (is.null(x$position) || length(x$position) == 0) {
    x$position <- 1:nLoc(x)
    if (verbose >= 2) {
      cat(warn("  Position slot is empty. Using a sequence from one to the number of loci in the dataset as dummy position.\n"))
    }
  }

  # Recode chromosome names to integer codes with a stable mapping:
  # codes are assigned in alphabetical order of the distinct names, so the
  # same name always gets the same code for a given name set, regardless of
  # factor level bookkeeping.
  chrom_names <- as.character(x$chromosome)
  chrom_levels <- sort(unique(chrom_names))
  x$chromosome <- as.factor(match(chrom_names, chrom_levels))
  if (verbose >= 2) {
    cat(report("  Chromosome name to code mapping:\n"))
    cat(report(paste0(
      "    ", chrom_levels, " -> ", seq_along(chrom_levels), "\n"
    )))
  }

  geno_tmp <- data.frame(rs = locNames(x),
                         alleles= x$loc.all,
                         chrom= x$chromosome,
                         pos= x$position,
                         strand="+",
                         assembly=NA,
                         center= NA,
                         protLSID= NA,
                         assayLSID= NA,
                         panel=NA,
                         QCcode=NA)

  res_output <- cbind(geno_tmp,xx)
  res_output <- as.matrix(res_output)
  res_output[] <- as.character(res_output)
  res_output <- as.matrix(rbind(colnames(res_output),res_output))
  res_output <- as.data.frame(res_output)

  # write the hapmap data.frame as tab-delimited text; the header row is
  # already the first data row, as GAPIT expects with header = FALSE
  utils::write.table(
    res_output,
    file = outfilespec,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  if (verbose >= 2) {
    cat(report("  Output file:", outfilespec, "\n"))
  }

  # FLAG SCRIPT END

  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }


  return(invisible(res_output))
}
