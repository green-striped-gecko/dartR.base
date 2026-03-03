#' @name gl2sa
#' @title Converts genlight objects to the format used in the SNPassoc package
#' @family linker

#' @description
#' This function exports a genlight object into a SNPassoc object and requires 
#' the package SNPassoc to be installed. 
#' 
#' To install this package check their github repository:
#' 
#' \url{https://github.com/isglobal-brge/SNPassoc} 
#' 
#' and/or their webpage:
#' 
#' \url{https://isglobal-brge.r-universe.dev/SNPassoc}
#' 
#'  and/or use:
#'  
#' \code{install.packages('SNPassoc',
#' repos = c('https://isglobal-brge.r-universe.dev',
#' 'https://cloud.r-project.org'))}
#' 
#' @param x Name of the genlight object containing the SNP data [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @author Bernd Guber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @references
#' Gonzalez, J.R., Armengol, L., Sol?, X., Guin?, E., Mercader, J.M., Estivill,
#' X. and Moreno, V. (2017). SNPassoc: an R package to perform whole genome
#' association studies. Bioinformatics 23:654-655.
#' @export
#' @return Returns an object of class 'snp' to be used with SNPassoc.

gl2sa <- function(x, 
                  verbose = NULL) {
  
  sa <- NULL
  
  # CHECK IF PACKAGES ARE INSTALLED
  if (!(requireNamespace("parallel", quietly = TRUE))) {
    stop(error(
      "Package parallel needed for this function to work. Please install it."
    ))
  }
  
  if (!(requireNamespace("pegas", quietly = TRUE))) {
    stop(error(
      "Package pegas needed for this function to work. Please install it."
    ))
  }
  
  if (!(requireNamespace("SNPassoc", quietly = TRUE))) {
    stop(
      error(
        "To use this function you need to install package: SNPassoc. Please refer to the help of the function for instructions (?gl2sa)."
      )
    )
  } else {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # DO THE JOB
    
    if (verbose >= 2) {
      cat(report("  Writing data to SNPassoc object\n"))
    }
    pop <- gl2gi(x)
    xxx <- pegas::as.loci(pop)[, -1]
    sa <- SNPassoc::setupSNP(data.frame(xxx), 1:ncol(xxx), )
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
      cat(report("Completed:", funname, "\n"))
    }
    
    return(sa)
  }
}
