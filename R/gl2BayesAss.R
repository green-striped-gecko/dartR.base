#' @name gl2bayesAss
#' @title Converts a genlight object into bayesAss (BA3) input format
#' @family linker
#' 
#' @description
#' This function exports a genlight object into bayesAss format and save it into a
#' file.
#' This function only caters for \code{ploidy=2}.
#' 
#' @param x Name of the genlight object containing the SNP data [required].
#' @param ploidy Set the ploidy [defaults 2].
#' @param outfile File name of the output file [default 'gl.BayesAss.txt'].
#' @param outpath Path where to save the output file [default global working 
#' directory or if not specified, tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' 
#' @author Custodian: Carlo Pacioni (Post to
#'  \url{https://groups.google.com/d/forum/dartr})
#'  
#' @examples
#' require("dartR.data")
#' #only the first 100 due to check time
#' gl2bayesAss(platypus.gl[,1:100], outpath=tempdir())
#' @references
#' Mussmann S. M., Douglas M. R., Chafin T. K. and Douglas M. E. (2019) BA3-SNPs: 
#' Contemporary migration reconfigured in BayesAss for next-generation sequence data. 
#' Methods in Ecology and Evolution 10, 1808-1813.
#' 
#' Wilson G. A. and Rannala B. (2003) Bayesian Inference of Recent Migration Rates 
#' Using Multilocus Genotypes. Genetics 163, 1177-1191.
#' 
#' @export
#' @return  returns the input file as data.table
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
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  if(ploidy != 2) stop("This function only caters for ploidy=2")
  
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
