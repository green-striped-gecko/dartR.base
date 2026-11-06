#' @name utils.ped2plink
#' @title Reads a PLINK .ped fileset with PLINK allele orientation
#'
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.
#'
#' @details
#' Reads a PLINK text fileset (.ped, optionally with a .map) in R with
#' \code{\link{utils.read.ped}}, with no PLINK executable required, and
#' returns the same three components as \code{snpStats::read.plink()} so
#' that \code{\link{gl.read.PLINK}} can treat text and binary input
#' identically.
#'
#' \code{utils.read.ped} codes alleles in order of first appearance in the
#' file. PLINK 1.x instead makes the minor allele A1 and the major allele
#' A2, and on a tie the first-encountered allele becomes A2. Alleles are
#' switched here per locus to match that rule, so the genotype dosage
#' (count of allele.2) is the same as that obtained from a .bed fileset
#' written by \code{plink --make-bed} from the same .ped.
#'
#' @param ped Path to the .ped file [required].
#' @param map Path to the matching .map file, or NULL when there is none.
#' Without a .map, loci are named locus.1, locus.2, ... and chromosome and
#' position are set to NA [default NULL].
#' @return A list with \code{genotypes} (a SnpMatrix), \code{fam} (a
#' data.frame with pedigree, member, father, mother, sex, affected) and
#' \code{map} (a data.frame with chromosome, snp.name, cM, position,
#' allele.1, allele.2, one row per locus, row names = snp.name).
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @keywords internal
#' @importFrom snpStats col.summary switch.alleles

utils.ped2plink <- function(ped, map = NULL) {
  if (is.null(map)) {
    x <- utils.read.ped(file = ped)
    m <- x$map
    map.df <- data.frame(
      chromosome = rep(NA_integer_, nrow(m)),
      snp.name = as.character(m$snp.name),
      cM = rep(NA_real_, nrow(m)),
      position = rep(NA_integer_, nrow(m)),
      allele.1 = as.character(m$allele.1),
      allele.2 = as.character(m$allele.2),
      stringsAsFactors = FALSE
    )
  } else {
    x <- utils.read.ped(file = ped, snps = map)
    m <- x$map
    # .map columns: chromosome, snp id, cM, bp
    map.df <- data.frame(
      chromosome = m[[1]],
      snp.name = as.character(m$snp.names),
      cM = m[[3]],
      position = m[[4]],
      allele.1 = as.character(m$allele.1),
      allele.2 = as.character(m$allele.2),
      stringsAsFactors = FALSE
    )
  }
  gen <- x$genotypes

  # Orient alleles as PLINK does: allele.2 = major allele; on a tie the
  # first-encountered allele (allele.1 in utils.read.ped) becomes allele.2.
  cs <- snpStats::col.summary(gen)
  n2 <- round(2 * cs$Calls * cs$RAF)
  n2[is.na(n2)] <- 0
  n1 <- 2 * cs$Calls - n2
  sw <- which(n1 >= n2)
  if (length(sw) > 0) {
    gen <- snpStats::switch.alleles(gen, sw)
    a1 <- map.df$allele.1
    map.df$allele.1[sw] <- map.df$allele.2[sw]
    map.df$allele.2[sw] <- a1[sw]
  }
  rownames(map.df) <- map.df$snp.name

  list(genotypes = gen, fam = x$fam, map = map.df)
}
