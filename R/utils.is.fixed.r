#' @name utils.is.fixed
#' @title An internal function to tests if two populations are fixed at a given locus
#' @family utilities
#' 
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.

#' @param s1 Percentage SNP allele or sequence tag frequency for the first 
#' population [required].
#' @param s2 Percentage SNP allele or sequence tag frequency for the second 
#' population [required].
#' @param tloc Threshold value for tolerance in when a difference is regarded as
#'  fixed [default 0].
#'  
#' @details
#' This script compares two percent allele frequencies
#' and reports TRUE if they represent a fixed difference, FALSE otherwise.
#' 
#'  A fixed difference at a locus occurs when two populations share no alleles,
#'  noting that SNPs are biallelic (ploidy=2).
#' Tolerance in the definition of a fixed difference is provided by the t
#'  parameter. For example, t=0.05 means that SNP allele frequencies of 95,5 and
#'  5,95 percent will be reported as fixed (TRUE).
#'  
#' @author Maintainer: Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})

#' @seealso \code{\link{gl.fixed.diff}}
#'  
#' @export
#' @return TRUE (fixed difference) or FALSE (alleles shared) or NA (one or both s1 or s2 missing)

# Examples for testing
# utils.is.fixed(s1=100, s2=0, tloc=0)
# utils.is.fixed(96, 4, tloc=0.05)

utils.is.fixed <- function(s1, 
                     s2, 
                     tloc = 0) {
    if (is.na(s1) | is.na(s2)) {
        result <- NA
    } else {
        # Transform to fixed, if outside tollerance
        if (s1 <= tloc * 100) {
            s1 <- 0
        }
        if (s1 >= 100 - tloc * 100) {
            s1 <- 100
        }
        if (s2 <= tloc * 100) {
            s2 <- 0
        }
        if (s2 >= 100 - tloc * 100) {
            s2 <- 100
        }
        # Score
        result <- 0
        if ((s1 == 0) & (s2 == 100)) {
            result <- 1
        }
        if ((s1 == 100) & (s2 == 0)) {
            result <- 1
        }
    }
    return(result)
}

# Test function utils.is.fixed(0,100) utils.is.fixed(100,0) utils.is.fixed(80,0) utils.is.fixed(100,NA) 
#utils.is.fixed(0,NA) utils.is.fixed(NA,0) utils.is.fixed(NaN,100)
# utils.is.fixed(NaN,0) utils.is.fixed(100,NaN) utils.is.fixed(0,NaN) utils.is.fixed(NaN,NaN) 
#utils.is.fixed(NA,NA) utils.is.fixed(50,50,tloc=0.05)
