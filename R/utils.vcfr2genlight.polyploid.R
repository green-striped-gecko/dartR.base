#' @name utils.vcfr2genlight.polyploid
#' @title Utility function to convert polyploid vcfR object as genlight
#' @family utilities
#' @description 
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.
#' @param x Name of the vcfR object [defined in function \code{\link{gl.read.vcf}}].
#' @param mode2 genotype: all heterozygous sites will be coded as 1 regardless ploidy level,
#' dosage: sites will be codes as copy number of alternate allele [default 'genotype'].
#' @param n.cores Number of cores [default 1]
#' @details
#' This function uses parameters from \code{\link{gl.read.vcf}} for conversion
#' The ploidy of each individual is derived from the allele count (arity) of
#' its genotype calls. Half-missing calls (e.g. 0/., ./1) are treated as
#' missing data (NA) in both modes.
#' Please carefully check the data if "dosage" mode is used. (codes were modified from
#' 'vcfR2genlight' in vcfR packge to convert polyploid data)
#' @author Custodian: Ching Ching Lau -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' datatype <- utils.vcfr2genlight.polyploid(x=vcfr, mode2="genotype")
#' }
#' @references
#' \itemize{
#' \item Knaus, B. J., & Grunwald, N. J. (2017). 
#' vcfr: a package to manipulate and visualize variant call format data in R. 
#' Molecular ecology resources, 17(1), 44-53.
#' \item Knaus, B. J., Grunwald, N. J., Anderson, E. C., Winter, D. J., 
#' Kamvar, Z. N., & Tabima, J. F. (2023). Package 'vcfR'.
#' \href{https://github.com/knausb/vcfR/blob/master/R/vcfR_conversion.R}{vcfR}
#' }
#' @return genlight object
#' @export

utils.vcfr2genlight.polyploid <- function(x, n.cores=1, mode2="genotype") {
    bi <- vcfR::is.biallelic(x)
    if(sum(!bi) > 0){
      msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
      msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
      msg <- c(msg, "\n", paste(sum(!bi), 'loci will be omitted from the genlight object.'))
      warning(msg)
      x <- x[bi,]
    }
    
    x <- vcfR::addID(x)
    
    CHROM <- x@fix[,'CHROM']
    POS   <- x@fix[,'POS']
    ID    <- x@fix[,'ID']
    
    x <- vcfR::extract.gt(x)

    # Ploidy is the allele count per GT call (arity), read off the raw
    # calls before the separators are stripped -- exact, unlike the
    # max-dosage inference adegenet applies at object creation (F2).
    ploidy_ind <- apply(x, 2, function(g) {
      g <- g[!is.na(g)]
      if (length(g) == 0) return(NA_integer_)
      max(nchar(gsub("[^/|]", "", g))) + 1L
    })
    if (all(is.na(ploidy_ind))) {
      ploidy_ind[] <- 2L
    } else if (any(is.na(ploidy_ind))) {
      # individuals with no called genotypes take the modal ploidy
      ploidy_ind[is.na(ploidy_ind)] <-
        as.integer(names(which.max(table(ploidy_ind))))
    }

    x <- gsub("/", "", x)
    x <- gsub("|", "", x, fixed = TRUE)
    # Half-missing calls (e.g. 0/., ./1) still carry a "." after the
    # separators are stripped; they are missing data, not evidence of a
    # heterozygote or a dosage -- convert them to NA before the mode
    # rules (F1).
    x[grepl(".", x, fixed = TRUE)] <- NA
    # code all polyploid heterozygous sites to 1
    if (mode2=="genotype"){
      x[stringr::str_count(as.character(x),"0") == nchar(as.character(x))] <- 0
      x[stringr::str_count(as.character(x),"1") == nchar(as.character(x))] <- 2
      x[nchar(as.character(x)) != 1 & stringr::str_count(as.character(x),"1")/nchar(as.character(x)) < 1] <- 1
    } else if (mode2=="dosage") {
      #allow different codes other than 0,1,2,NA
      # all 0
      x[stringr::str_count(as.character(x),"0") == nchar(as.character(x))] <- 0
      # all 1
      x[which(stringr::str_count(as.character(x),"1") == nchar(as.character(x)))] <- 
        nchar(x[which(stringr::str_count(as.character(x),"1") == nchar(as.character(x)))])
      # heterozygous
      x[which(nchar(as.character(x)) != 1 & stringr::str_count(as.character(x),"1")/nchar(as.character(x)) < 1 &
                stringr::str_count(as.character(x),"1")/nchar(as.character(x)) > 0)] <-
        stringr::str_count(x[which(nchar(as.character(x)) != 1 & 
                                     stringr::str_count(as.character(x),"1")/nchar(as.character(x)) < 1 &
                                     stringr::str_count(as.character(x),"1")/nchar(as.character(x)) > 0)],"1")
    } else {
      cat(error("  Please choose 'genotype' or 'dosage' mode \n"))
      stop()
    }
    #  dim(x)
    if( requireNamespace('adegenet') ){
      x <- new('genlight', t(x), n.cores=n.cores)
    } else {
      warning("adegenet not installed")
    }
    #  x <- adegenet::as.genlight(t(x), n.cores=3)
    #  x <- adegenet::as.genlight(t(x))
    # Stamp the arity-derived ploidy over adegenet's max-dosage
    # inference (F2).
    adegenet::ploidy(x) <- ploidy_ind
    adegenet::chromosome(x) <- CHROM
    adegenet::position(x)   <- POS
    adegenet::locNames(x)   <- ID
    
    return(x)
  }