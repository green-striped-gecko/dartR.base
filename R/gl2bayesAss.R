#' @name gl2bayesAss
#' @title Converts a genlight object into bayesAss (BA3) input format
#' @family linker
#'
#' @description
#' This function exports a genlight object into bayesAss format and saves it
#' into a file.
#' This function only caters for \code{ploidy = 2}.
#'
#' @details
#' The output file has one row per individual per locus, in the format
#' expected by BayesAss (BA3) and BA3-SNPs: individual identifier,
#' population, locus name, then the two alleles coded as 1 (reference) and
#' 2 (alternate). Missing genotypes are written as 0 0.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param ploidy Confirmation that the data are diploid; 2 is the only
#' accepted value, other values stop the function [default 2].
#' @param outfile File name of the output file [default 'gl.BayesAss.txt'].
#' @param outpath Path where to save the output file [default global working
#' directory or if not specified, tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @return returns the input file as data.table
#'
#' @author Author(s): Carlo Pacioni. Custodian: Carlo Pacioni -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @references
#' Mussmann S. M., Douglas M. R., Chafin T. K. and Douglas M. E. (2019) BA3-SNPs:
#' Contemporary migration reconfigured in BayesAss for next-generation sequence data.
#' Methods in Ecology and Evolution 10, 1808-1813.
#'
#' Wilson G. A. and Rannala B. (2003) Bayesian Inference of Recent Migration Rates
#' Using Multilocus Genotypes. Genetics 163, 1177-1191.
#'
#' @examples
#' require("dartR.data")
#' #only the first 100 due to check time
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' gl2bayesAss(platypus.gl[,1:100], outpath=tempdir())
#' @export
gl2bayesAss <-  function(x, 
                         ploidy=2, 
                         outfile="gl.BayesAss.txt", 
                         outpath=NULL,
                         verbose=NULL) {
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
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (ploidy != 2) {
    stop(error("Fatal Error: this function only caters for ploidy = 2\n"))
  }
  
  # DO THE JOB
  # Set NULL to variables to pass CRAN checks
  Locus<-Pop<-i.V2<-i.V3<-locus<-NULL
  
  m <- as.matrix(x)
  dtref<- data.table(c(0,1,2,NA), c(1,1,2,0), c(1,2,2,0))
  mdt <- data.table(m, keep.rownames = T)
  
  byLoc <- function(i, mdt, dtref) {
    loc <- names(mdt)[1 + i]
    mdtsub <- mdt[,c("rn", loc), with=F]
    setnames(mdtsub, loc, "locus")
    mdtsub[dtref, on=c(locus="V1"), c("All1", "All2") := list(i.V2, i.V3)]
    mdtsub[, locus := NULL]
    mdtsub[, Locus := loc]
    mdtsub[, Pop := pop(x)]
    setcolorder(mdtsub, c("rn", "Pop", "Locus", "All1", "All2"))
    return(mdtsub)
  }
  
  l <- lapply(seq_len(ncol(m)), byLoc, mdt, dtref)
  res <- rbindlist(l)
  setkeyv(res, cols=c("Pop", "rn"))
  write.table(res, file.path(outpath, outfile), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  
  return(invisible(res))
  }

