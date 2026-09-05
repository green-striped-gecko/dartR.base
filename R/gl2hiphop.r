#' @name gl2hiphop
#' @title Converts a genlight objects into hiphop format
#' @family linker

#' @description
#' This function exports genlight objects to the format used by the parentage
#' assignment R package hiphop. Hiphop can be used for paternity and maternity
#' assignment and outperforms conventional methods where closely related
#' individuals occur in the pool of possible parents. The method compares the
#' genotypes of offspring with any combination of potentials parents and scores
#' the number of mismatches of these individuals at bi-allelic genetic markers
#'  (e.g. Single Nucleotide Polymorphisms).

#' @details
#' The dartR and hiphop genotype encodings differ: dartR scores
#' 0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate,
#' while hiphop expects 0 = homozygote, 1 = the other homozygote,
#' 2 = heterozygote. The conversion therefore maps 0 -> 0, 1 -> 2,
#' 2 -> 1, and missing values (NA) are preserved.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#'
#' @return Dataframe containing all the genotyped individuals (offspring and
#'  potential parents) and their genotypes scored using bi-allelic markers.
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' result <- gl2hiphop(testset.gl)
#' }
#'
#' @references
#' Cockburn, A., Penalba, J.V.,Jaccoud, D.,Kilian, A., Brouwer, L.,Double, M.C.,
#'  Margraf, N., Osmond, H.L., van de Pol, M. and Kruuk, L.E.B.(in revision).
#'  HIPHOP: improved paternity assignment among close relatives using a simple
#'  exclusion method for bi-allelic markers. Molecular Ecology Resources, DOI to
#'  be added upon acceptance
#'
# The magrittr pipe import is retained here although this function no
# longer uses it: gl.map.interactive, gl.pcoa.plot and utils.heatmap use
# %>% at run time and this was the only @importFrom declaring it.
# Follow-up: move the declaration to a package-level file.
#' @importFrom dplyr %>%
#' @export

gl2hiphop <- function(x,
                      verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)
    
    # CHECK DATATYPE
    # hiphop scoring is defined for biallelic SNP dosages only (F1)
    datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

    # DO THE JOB

    # swap the 1 and 2 codes numerically: dartR 0/1/2 (hom ref/het/hom
    # alt) -> hiphop 0/2/1 (hom/het/other hom); NA is preserved (F4)
    m <- as.matrix(x)
    m[] <- c(0, 2, 1)[m + 1]
    x <- as.data.frame(m)

    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x)
}
