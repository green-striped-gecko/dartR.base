# Estimates the lower bound of the number of undetected private alleles
# between a pair of populations (Good-Turing / Chao, equation 2c of Chao
# et al. 2017), separately for each direction of the comparison.
#
# The estimate is based on the pair only: for each locus carrying a
# private allele, the number of copies of that allele observed in the
# pooled sample of the two populations is counted directly, and f1/f2
# are the numbers of private alleles observed exactly once/twice.
# SNP data only (the caller reports NA for SilicoDArT).

utils.pa.Chao <- function(pop1_m, pop2_m) {
  p1_m <- as.matrix(pop1_m)
  p2_m <- as.matrix(pop2_m)
  p1alf <- colMeans(p1_m, na.rm = TRUE) / 2
  p2alf <- colMeans(p2_m, na.rm = TRUE) / 2
  pooled <- rbind(p1_m, p2_m)

  chao.dir <- function(host_alf, other_alf) {
    # loci where the host population carries an allele absent in the other
    alt.private <- which(other_alf == 0 & host_alf != 0)
    ref.private <- which(other_alf == 1 & host_alf != 1)
    n <- length(alt.private) + length(ref.private)
    if (n == 0) {
      return(0)
    }
    # copies of the private allele observed in the pooled pair sample
    alt.count <- colSums(pooled[, alt.private, drop = FALSE], na.rm = TRUE)
    ref.nonNA <- colSums(!is.na(pooled[, ref.private, drop = FALSE]))
    ref.count <- 2 * ref.nonNA -
      colSums(pooled[, ref.private, drop = FALSE], na.rm = TRUE)
    counts <- c(alt.count, ref.count)
    f1 <- sum(counts == 1)
    f2 <- sum(counts == 2)
    if (f2 == 0) {
      # bias-corrected form when no doubletons are observed
      ((n - 1) / n) * f1 * (f1 - 1) / 2
    } else {
      ((n - 1) / n) * (f1 ^ 2) / (2 * f2)
    }
  }

  return(list(chao.dir(p1alf, p2alf), chao.dir(p2alf, p1alf)))
}
