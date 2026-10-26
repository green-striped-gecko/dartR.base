#' @name utils.vcf.depth
#' @title Utility function to read locus depth and allele read counts from a
#' vcfR object
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY
#' END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.
#' @param vcf A vcfR object [defined in function \code{\link{gl.read.vcf}}].
#' @param loci Indices of the records kept as loci, in the order of the loci
#' in the genlight object [default all records].
#' @param read.counts If TRUE, the reference and alternative read counts of
#' each genotype are also returned [default FALSE].
#' @param chunk Number of records processed at a time, which bounds the
#' memory used [default 10000].
#' @details
#' rdepth is the mean of FORMAT DP over the called genotypes of a locus.
#' Missing and half-missing calls (e.g. ./., 0/.) are not called, as in
#' \code{\link{utils.vcfr2genlight.polyploid}}; some callers report DP = 0
#' for them. A locus without called genotypes or without DP gets NA.
#'
#' Read counts come from FORMAT AD (first value REF, second ALT). Where a
#' genotype has no AD, they come from FORMAT BSDP (read counts of the bases A,
#' C, G and T, as written by some callers) taken at the REF and ALT bases;
#' records whose REF or ALT is not a single base get NA. Counts are stored as
#' reported, also for uncalled genotypes.
#' @return A list with rdepth (one value per locus, NULL when the vcf has no
#' DP) and, when read.counts is TRUE, ref.count and alt.count (integer
#' matrices, individuals x loci, NULL when the vcf has neither AD nor BSDP).
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' depth <- utils.vcf.depth(vcf, read.counts = TRUE)
#' }
#' @export

utils.vcf.depth <- function(vcf,
                            loci = seq_len(nrow(vcf@fix)),
                            read.counts = FALSE,
                            chunk = 10000L) {
  res <- list(rdepth = NULL, ref.count = NULL, alt.count = NULL)
  if (ncol(vcf@gt) < 2 || length(loci) == 0) {
    return(res)
  }

  fields <- unique(unlist(strsplit(unique(vcf@gt[loci, "FORMAT"]), ":",
                                   fixed = TRUE)))
  has.dp <- "DP" %in% fields
  has.ad <- read.counts && "AD" %in% fields
  has.bsdp <- read.counts && "BSDP" %in% fields
  n.ind <- ncol(vcf@gt) - 1L
  n.loc <- length(loci)

  if (has.dp) {
    res$rdepth <- rep(NA_real_, n.loc)
  }
  if (has.ad || has.bsdp) {
    res$ref.count <- matrix(NA_integer_, n.ind, n.loc)
    res$alt.count <- matrix(NA_integer_, n.ind, n.loc)
  }
  if (!has.dp && is.null(res$ref.count)) {
    return(res)
  }

  bases <- c("A", "C", "G", "T")
  for (s in seq(1L, n.loc, by = chunk)) {
    cols <- s:min(n.loc, s + chunk - 1L)
    v <- vcf[loci[cols], ]

    if (has.dp) {
      gt <- vcfR::extract.gt(v, element = "GT", IDtoRowNames = FALSE)
      dp <- vcfR::extract.gt(v, element = "DP", as.numeric = TRUE,
                             IDtoRowNames = FALSE)
      dp[is.na(gt) | grepl(".", gt, fixed = TRUE)] <- NA
      n.called <- rowSums(!is.na(dp))
      res$rdepth[cols] <- ifelse(n.called > 0,
                                 rowSums(dp, na.rm = TRUE) / n.called, NA)
    }

    if (!is.null(res$ref.count)) {
      ref <- alt <- matrix(NA_real_, length(cols), n.ind)
      if (has.ad) {
        ad <- vcfR::extract.gt(v, element = "AD", IDtoRowNames = FALSE)
        ref <- vcfR::masplit(ad, record = 1, sort = 0)
        alt <- vcfR::masplit(ad, record = 2, sort = 0)
      }
      # genotypes without AD take the counts of their REF and ALT bases
      need <- is.na(ref) & is.na(alt)
      if (has.bsdp && any(need)) {
        bsdp <- vcfR::extract.gt(v, element = "BSDP", IDtoRowNames = FALSE)
        ref.base <- match(v@fix[, "REF"], bases)
        alt.base <- match(v@fix[, "ALT"], bases)
        # at an indel the base counts do not separate the two alleles, even
        # when one allele is a single base
        snv <- !is.na(ref.base) & !is.na(alt.base)
        ref.base[!snv] <- NA
        alt.base[!snv] <- NA
        for (k in seq_along(bases)) {
          if (!any(ref.base %in% k | alt.base %in% k)) next
          count <- vcfR::masplit(bsdp, record = k, sort = 0)
          # a vector of one value per record recycles down the columns
          sel <- need & (ref.base %in% k)
          ref[sel] <- count[sel]
          sel <- need & (alt.base %in% k)
          alt[sel] <- count[sel]
        }
      }
      res$ref.count[, cols] <- as.integer(t(ref))
      res$alt.count[, cols] <- as.integer(t(alt))
    }
  }

  res
}
