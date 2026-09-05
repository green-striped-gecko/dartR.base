#' @name gl.read.fasta
#' @title Reads FASTA files and converts them to genlight object
#' @family io

#' @description
#' The following IUPAC Ambiguity Codes are taken as heterozygotes:
#' \itemize{
#'  \item M is heterozygote for	AC and CA	
#'  \item R is heterozygote	for AG and GA	
#'  \item W is heterozygote	for AT and TA	
#'  \item S is heterozygote	for CG and GC
#'  \item Y is heterozygote	for CT and TC
#'  \item K is heterozygote	for GT and TG
#'  }

#' The following IUPAC Ambiguity Codes are taken as missing data:
#' \itemize{
#'  \item V
#'  \item H 
#'  \item D 
#'  \item B
#'  \item N 
#'  }

#'  The function can deal with missing data in individuals, e.g. when FASTA 
#'  files have different number of individuals due to missing data.

#'  The allele with the highest frequency is taken as the reference allele.

#'  SNPs with more than two alleles are skipped.

#' @param fasta.files Fasta files to read [required].
#' @param parallel A logical indicating whether multiple cores -if available-
#'  should be used for the computations (TRUE), or not (FALSE); requires the
#'   package parallel to be installed [default FALSE].
#' @param n.cores If parallel is TRUE, the number of cores to be used in the
#'  computations; if NULL, then the maximum number of cores available on the
#'   computer is used [default NULL].
#' @param fbm Logical. If TRUE, the returned genlight object will contain a
#' file-backed matrix (fbm) in its @genome slot. This is useful for very
#' large datasets that do not fit into RAM. Note that using fbm objects
#' requires the package bigsnpr to be installed. [default FALSE]. to back convert 
#' use \code{gl.fbm2gen()}.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' 
#' @details
#' Ambiguity characters are often used to code heterozygotes. However, using
#'  heterozygotes as ambiguity characters may bias many estimates. See more
#'   information in the link below:
#' \url{https://evodify.com/heterozygotes-ambiguity-characters/}
#'
#' Each record in a FASTA file must occupy exactly two lines: a '>' header
#' line followed by the full sequence on a single line. Line-wrapped
#' (multi-line) FASTA is not supported and is rejected with an error.
#' 
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' # Folder where the fasta files are located. 
#' folder_samples <- system.file('extdata', package ='dartR.data')
#' # listing the FASTA files, including their path. Files have an extension
#' # that contains "fas".
#' file_names <- list.files(path = folder_samples, pattern = "*.fas", 
#' full.names = TRUE)
#' # reading fasta files
#' obj <- gl.read.fasta(file_names)
#'   
#' @export
#' @return A genlight object.

gl.read.fasta <- function(fasta.files,
                          parallel = FALSE,
                          n.cores = NULL,
                          fbm=FALSE,
                          verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING

  # The reader requires each record to occupy exactly two lines (a '>'
  # header and the full sequence on a single line); line-wrapped FASTA
  # mis-groups records silently, so reject it before reading
  for (f in fasta.files) {
    tmp <- scan(f, what = "character", sep = "\n", quiet = TRUE)
    headers <- grep("^>", tmp)
    if (length(headers) == 0 ||
        length(tmp) != 2 * length(headers) ||
        !all(headers == seq(1, by = 2, length.out = length(headers)))) {
      stop(error(
        "Fatal Error: file", basename(f), "is not two-line FASTA;",
        "each record must be a '>' header line followed by the full sequence on a single line\n"
      ))
    }
  }

  # DO THE JOB

  gl_list <- lapply(fasta.files,
                    utils.read.fasta,
                    parallel = parallel,
                    n.cores = n.cores,
                    verbose = verbose)

  # A file without polymorphism yields NULL from utils.read.fasta; if no
  # file yielded SNPs there is nothing to build an object from
  if (all(vapply(gl_list, is.null, logical(1)))) {
    stop(error(
      "Fatal Error: no SNPs found; the input file(s) contain no polymorphic positions\n"
    ))
  }

  x <- merge_gl_fasta(gl_list, parallel = parallel, verbose = verbose)
  
  x <- gl.compliance.check(x, verbose = verbose)
 
   # add history
  x@other$history <- list(match.call())
  x <- gl.recalc.metrics(x, verbose = 0)
  
  if (verbose >= 2) {
    cat(
      important(
        "  Genlight object does not have individual metrics. You need to add them 'manually' to the @other$ind.metrics slot.\n"
      )
    )
  }
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  #convert to fbm, only when requested
  if (fbm) {
    x <- gl.gen2fbm(x, verbose = verbose)
    if (verbose > 2) {
      cat(report(" Created a file-backed matrix (fbm) dartR object\n"))
    }
  } else {
    x@fbm <- NULL
  }
  
  # RETURN
  return(invisible(x))
  
}
