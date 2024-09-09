#' @name gl2geno
#' @title Converts a genlight object to geno format from package LEA
#' @family linker

#' @description
#' The function converts a genlight object (SNP or presence/absence
#'  i.e. SilicoDArT data) into a file in the 'geno' and the 'lfmm' formats from 
#'  (package LEA).
#'  
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file [default 'gl_gapit'].
#' @param outpath Path where to save the output file [default global working 
#' directory or if not specified, tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' 
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#' t1 <- platypus.gl
#' # assign chromosome
#' t1$chromosome <- t1$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1
#' # assign position
#' t1$position <- t1$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
#' res <- gl2gapit(t1)
#' 
#' @export
#' @return  returns no value (i.e. NULL)

gl2gapit <- function(x, 
                     outfile = "gl_gapit",
                     outpath = NULL,
                     verbose = NULL){
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # SET WORKING DIRECTORY
  outpath <- gl.check.wd(outpath,verbose=0)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # DO THE JOB
  
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
  xx <- gsub("/", "", xx)
  xx <- as.data.frame(xx)
  xx <- t(xx)
  colnames(xx) <- indNames(x)
  
  x$chromosome <- as.factor(as.numeric(x$chromosome))
  
  geno_tmp <- data.frame(rs = locNames(x),
                         alleles= x$loc.all,
                         chrom= x$chromosome,
                         pos= x$position,
                         strand="+",
                         assembly="Oilpalm",
                         center= NA,
                         protLSID= NA,
                         assayLSID= NA,
                         panel=NA,
                         QCcode=NA)
  
  res_output <- cbind(geno_tmp,xx)
  res_output <- as.matrix(res_output)
  res_output[] <- as.character(res_output)
  res_output <- as.matrix(rbind(colnames(res_output),res_output))
  
  
  if (verbose > 0) {
    cat(report("  Output files:", paste(
      paste0(outfile, ".geno", ".lfmm."), sep = ""
    ), "\n"))
  }
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  

  return(invisible(res_output))
}
