#' @name utils.impute
#' @title An internal script [Custodian to provide a title]
#' @family utilities
#' 
#' @description 
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.

#' @param snp_matrix [Custodian to provide parameter description]
#' @param parallel [Custodian to provide parameter description]

#' @details
#' #[Custodian to provide details for future you]
#' 
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})

# @export
#' @return The resultant genlight object


matrix2gen <- function(snp_matrix, parallel = FALSE) {
  if (parallel) {
    i@gen <-
      parallel::mclapply(1:nrow(snp_matrix), function(i)
        new("SNPbin", as.integer(snp_matrix[i, ])), mc.silent = TRUE, mc.cleanup =
          TRUE, mc.preschedule = FALSE)
  } else {
    lapply(1:nrow(snp_matrix), function(i)
      new("SNPbin", as.integer(snp_matrix[i, ])))
  }
}

# function to sample alleles using allele frequencies as probability
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

# function to sample genotypes based on Hardy-Weinberg equation
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

# Stream-impute an FBM.code256 column-block-wise, per population.
# Used by gl.impute() to avoid materialising the full nInd x nLoc matrix.
# Returns list(all_nas, imputations) as named integer vectors per population.
.stream_impute_pop_fbm <- function(fbm, ind_by_pop, method, chunk = 2048L) {
  stopifnot(method %in% c("frequency", "HW", "random"))
  n_loc <- ncol(fbm)
  starts <- seq.int(1L, n_loc, by = chunk)
  all_nas     <- setNames(integer(length(ind_by_pop)), names(ind_by_pop))
  imputations <- setNames(integer(length(ind_by_pop)), names(ind_by_pop))

  for (s in starts) {
    e <- min(s + chunk - 1L, n_loc)
    block <- fbm[, s:e, drop = FALSE]

    for (pn in names(ind_by_pop)) {
      ri <- ind_by_pop[[pn]]
      if (!length(ri)) next
      sub <- block[ri, , drop = FALSE]

      col_all_na <- colSums(is.na(sub)) == length(ri)
      all_nas[[pn]] <- all_nas[[pn]] + sum(col_all_na)

      loc_na <- which(is.na(sub), arr.ind = TRUE)
      if (!nrow(loc_na)) next

      if (method == "random") {
        sub[loc_na] <- sample(0:2, size = nrow(loc_na), replace = TRUE)
        imputations[[pn]] <- imputations[[pn]] + nrow(loc_na)
      } else {
        q <- colMeans(sub, na.rm = TRUE) / 2
        qv <- q[loc_na[, 2]]
        ok <- !is.na(qv)
        if (any(ok)) {
          fill <- if (method == "frequency") {
            as.integer(round(2 * qv[ok]))
          } else {
            vapply(qv[ok], function(qf) as.integer(sample_genotype(q_freq = qf)), integer(1))
          }
          sub[loc_na[ok, , drop = FALSE]] <- fill
          imputations[[pn]] <- imputations[[pn]] + sum(ok)
        }
      }
      block[ri, ] <- sub
    }
    # FBM.code256 stores NA as raw byte 3; writing NA_real_ via [<- silently
    # corrupts it to 0 (nan -> 0 coercion warning).
    block[is.na(block)] <- 3
    fbm[, s:e] <- block
  }
  list(all_nas = all_nas, imputations = imputations)
}

# Stream-impute FBM.code256 via nearest neighbour. D is the precomputed
# nInd x nInd distance matrix (with Inf on the diagonal and for undefined
# entries). Returns list(unfilled <per-ind>, all_nas_per_pop <named>).
.stream_impute_neighbour_fbm <- function(fbm, D, ind_by_pop, chunk = 2048L) {
  n_ind <- nrow(fbm)
  n_loc <- ncol(fbm)
  starts <- seq.int(1L, n_loc, by = chunk)
  ord_list <- lapply(seq_len(n_ind), function(i) order(D[i, ]))
  unfilled_per_ind <- integer(n_ind)
  all_nas_per_pop  <- setNames(integer(length(ind_by_pop)), names(ind_by_pop))

  for (s in starts) {
    e <- min(s + chunk - 1L, n_loc)
    block <- fbm[, s:e, drop = FALSE]

    for (pn in names(ind_by_pop)) {
      ri <- ind_by_pop[[pn]]
      if (!length(ri)) next
      sub <- block[ri, , drop = FALSE]
      all_nas_per_pop[[pn]] <- all_nas_per_pop[[pn]] +
        sum(colSums(is.na(sub)) == length(ri))
    }

    for (i in seq_len(n_ind)) {
      miss <- is.na(block[i, ])
      if (!any(miss)) next
      for (j in ord_list[[i]]) {
        if (!any(miss)) break
        cand <- block[j, miss]
        ok <- !is.na(cand)
        if (any(ok)) {
          miss_pos <- which(miss)
          fill_pos <- miss_pos[ok]
          block[i, fill_pos] <- cand[ok]
          miss[fill_pos] <- FALSE
        }
      }
      unfilled_per_ind[i] <- unfilled_per_ind[i] + sum(miss)
    }
    # FBM.code256 stores NA as raw byte 3; writing NA_real_ via [<- silently
    # corrupts it to 0 (nan -> 0 coercion warning).
    block[is.na(block)] <- 3
    fbm[, s:e] <- block
  }
  list(unfilled = unfilled_per_ind, all_nas_per_pop = all_nas_per_pop)
}

# Stream residual-NA fill on an FBM using stochastic s_alleles() draws
# parameterised by per-locus global allele frequencies (computed per block).
.stream_fill_residual_fbm <- function(fbm, chunk = 2048L) {
  n_loc <- ncol(fbm)
  starts <- seq.int(1L, n_loc, by = chunk)
  n_filled <- 0L
  for (s in starts) {
    e <- min(s + chunk - 1L, n_loc)
    block <- fbm[, s:e, drop = FALSE]
    loc_na <- which(is.na(block), arr.ind = TRUE)
    if (!nrow(loc_na)) next
    q <- colMeans(block, na.rm = TRUE) / 2
    qv <- q[loc_na[, 2]]
    ok <- !is.na(qv)
    if (any(ok)) {
      fill <- vapply(qv[ok], function(qf) as.integer(s_alleles(q_freq = qf)), integer(1))
      block[loc_na[ok, , drop = FALSE]] <- fill
      n_filled <- n_filled + sum(ok)
    }
    # FBM.code256 stores NA as raw byte 3; writing NA_real_ via [<- silently
    # corrupts it to 0 (nan -> 0 coercion warning).
    block[is.na(block)] <- 3
    fbm[, s:e] <- block
  }
  n_filled
}