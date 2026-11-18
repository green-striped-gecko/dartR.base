# Estimates the lower bound of the number of undetected private alleles
# between a pair of populations (Good-Turing / Chao, equation 2c of Chao
# et al. 2017), separately for each direction of the comparison.
#
# The estimate is based on the pair only: for each locus carrying a
# private allele, the number of copies of that allele observed in the
# pooled sample of the two populations is counted directly, and f1/f2
# are the numbers of private alleles observed exactly once/twice. The
# sample size n in the (n - 1) / n factor is the number of allele copies
# sampled in the pooled pair (twice the mean number of genotyped
# individuals at the private loci), not the number of private alleles.
# SNP data only (the caller reports NA for SilicoDArT).
#
# pop1_m, pop2_m: the two populations (genlight objects or genotype
# matrices, 0/1/2 dosages, individuals in rows).

utils.pa.Chao <- function(pop1_m, pop2_m) {
  p1_m <- as.matrix(pop1_m)
  p2_m <- as.matrix(pop2_m)
  .pa_chao(p1alf = colMeans(p1_m, na.rm = TRUE) / 2,
           p2alf = colMeans(p2_m, na.rm = TRUE) / 2,
           n.pooled = colSums(!is.na(p1_m)) + colSums(!is.na(p2_m)),
           s.pooled = colSums(p1_m, na.rm = TRUE) +
             colSums(p2_m, na.rm = TRUE))
}

# The estimates of utils.pa.Chao from per-locus terms of the pair, so that
# gl.report.pa can build them in blocks of loci without copying the
# populations: allele frequencies of each population (p1alf, p2alf), and,
# for the pooled pair, the number of genotyped individuals (n.pooled) and
# the sum of the 0/1/2 dosages (s.pooled). These sums are whole numbers, so
# adding the two populations' sums equals summing the pooled genotypes.
.pa_chao <- function(p1alf, p2alf, n.pooled, s.pooled) {
  chao.dir <- function(host_alf, other_alf) {
    # loci where the host population carries an allele absent in the other
    alt.private <- which(other_alf == 0 & host_alf != 0)
    ref.private <- which(other_alf == 1 & host_alf != 1)
    if (length(alt.private) + length(ref.private) == 0) {
      return(0)
    }
    # sample size: allele copies sampled in the pooled pair, averaged over
    # the private loci (loci differ in missing data)
    n <- mean(2 * n.pooled[c(alt.private, ref.private)])
    # copies of the private allele observed in the pooled pair sample
    alt.count <- s.pooled[alt.private]
    ref.count <- 2 * n.pooled[ref.private] - s.pooled[ref.private]
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
