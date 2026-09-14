# Internal helpers used by gl.impute, utils.transpose and the readers.
# Custodian: Luis Mijangos (Post to
# https://groups.google.com/d/forum/dartr)

#' Converts a genotype matrix (individuals x loci, 0/1/2/NA) to a list of
#' SNPbin objects suitable for the gen slot of a genlight object
#' @param snp_matrix Genotype matrix, one row per individual
#' @param parallel If TRUE, use parallel::mclapply (serial on Windows)
#' @return A list of SNPbin objects, one per row
#' @noRd
matrix2gen <- function(snp_matrix, parallel = FALSE) {
  if (parallel) {
    parallel::mclapply(1:nrow(snp_matrix), function(i)
      new("SNPbin", as.integer(snp_matrix[i, ])), mc.silent = TRUE,
      mc.cleanup = TRUE, mc.preschedule = FALSE)
  } else {
    lapply(1:nrow(snp_matrix), function(i)
      new("SNPbin", as.integer(snp_matrix[i, ])))
  }
}

#' Samples a genotype (0/1/2) by drawing two alleles with the alternate
#' allele frequency q_freq; NA in, NA out
#' @noRd
s_alleles <- function(q_freq) {
  if (is.na(q_freq)) {
    return(NA)
  }
  alleles_sampled <-
    paste0(sample(
      c("a", "A"),
      size = 2,
      prob = c(q_freq, 1 - q_freq),
      replace = T
    ), collapse = "")
  
  if (alleles_sampled == "AA") {
    alleles_sam <- 0
  }
  
  if (alleles_sampled == "aA") {
    alleles_sam <- 1
  }
  
  if (alleles_sampled == "Aa") {
    alleles_sam <- 1
  }
  
  if (alleles_sampled == "aa") {
    alleles_sam <- 2
  }
  
  return(as.numeric(alleles_sam))
}

#' Samples a genotype (0/1/2) from Hardy-Weinberg proportions at the
#' alternate allele frequency q_freq; NA in, NA out
#' @noRd
sample_genotype <- function(genotype_list = c(0, 1, 2), q_freq) {
  if (is.na(q_freq)) {
    return(NA)
  }
  #genotype probabilities based on Hardy-Weinberg equation
  # p^2 + 2pq + q^2 = 1
  geno_probs <- c(((1 - q_freq) ^ 2), # homozygote for the reference allele
                  (2 * (1 - q_freq) * q_freq), # heterozygote
                  (q_freq ^ 2)) # homozygote for the alternative allele)
  genotype_sampled <-
    sample(genotype_list, size = 1, prob = geno_probs)
  return(genotype_sampled)
}