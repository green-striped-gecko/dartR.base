# Het calculations
# df: genotype matrix, individuals x loci. ploidy: one value per individual
# (row), or a single value for all; the default 2 is the diploid calculation
# used by gl.report.heterozygosity.
pop.het_fun <- function(df,
                        n.invariant,
                        aHet,
                        bootstrap=TRUE,
                        ploidy = 2) {
  # rm loci that are all NA
  # otherwise these loci get Ho=0 which would not be correct
  loc.allNA <- colSums(is.na(df)) == nrow(df)
  df <- df[, !loc.allNA, drop = FALSE]
  
  diploid <- all(ploidy == 2)
  if (diploid) {
    Ho.loc <- colMeans(df == 1, na.rm = TRUE)
    q_freq <- colMeans(df, na.rm = TRUE) / 2
  } else {
    # Dosage data (0..k copies of the alternative allele). Ho is the gametic
    # heterozygosity, the probability that two allele copies drawn without
    # replacement from one individual differ, d(k - d) / choose(k, 2)
    # (Moody et al. 1993); for k = 2 it is the heterozygote indicator. q is
    # the alternative-allele share of all sampled allele copies.
    k <- rep_len(ploidy, nrow(df))
    Ho.loc <- colMeans(df * (k - df) / choose(k, 2), na.rm = TRUE)
    q_freq <- colSums(df, na.rm = TRUE) / colSums((!is.na(df)) * k)
  }
  n_loc.sample <- apply(df, 1, function(y) {
    sum(!is.na(y))
  })
  n_loc <- ncol(df)
  p_freq <- 1 - q_freq
  He.loc <- 2 * p_freq * q_freq
  n_ind.loc <- apply(df, 2, function(y) {
    sum(!is.na(y))
  })
  ### CP ### Unbiased He (i.e. corrected for sample size): N / (N - 1), with
  # N the number of sampled allele copies at the locus (2n for diploids)
  if (diploid) {
    n_all.loc <- 2 * as.numeric(n_ind.loc)
  } else {
    n_all.loc <- colSums((!is.na(df)) * k)
  }
  uHe.loc <- (n_all.loc / (n_all.loc - 1)) * He.loc
  
  FIS.loc <- 1 - (Ho.loc / uHe.loc)
  
  if(aHet) {
    all.res <- c(
      Ho.adj = mean(Ho.loc) * n_loc / (n_loc + n.invariant),
      Hexp.adj = mean(He.loc) * n_loc / (n_loc + n.invariant)
    )
  } else {
    all.res <- c(
      Ho = mean(Ho.loc, na.rm = TRUE),
      He = mean(He.loc, na.rm = TRUE),
      uHe = mean(uHe.loc, na.rm = TRUE),
      FIS = mean(FIS.loc, na.rm = TRUE)
    )
  }
  
  if(bootstrap) {
    return(all.res)
  } else {
    list(means=all.res, 
         byloc=list(Ho.loc=Ho.loc, He.loc=He.loc, uHe.loc=uHe.loc, FIS.loc=FIS.loc))
  }
  
}

compute.variability <- function(all.het, what.st=c("sd", "std.error"), what.het) {
  st.type <- match.arg(what.st)
  sapply(all.het, function(pop, st=st.type, het.type=what.het) {
    switch(st,
        sd=sd(pop[["byloc"]][[what.het]], na.rm = TRUE),
        std.error=std.error(pop[["byloc"]][[what.het]])
    )
  })
}

# bootstrapping function
pop.het <- function(df,
                    indices,
                    n.invariant = 0,
                    boot_method = "loc",
                    aHet=FALSE,
                    ploidy = 2) {
  
  df <- df[indices, , drop = FALSE]
  
  if(boot_method == "loc"){
    df <- t(df)
  }

  res <- pop.het_fun(df,
                     n.invariant = n.invariant,
                     aHet = aHet,
                     ploidy = ploidy)
  
  return(res)
  
}

# Counting individuals function
ind.count <- function(x) {
  # the loci that are completely missing
  loci.na <-
    which(colSums(is.na(as.matrix(x))) == nrow(as.matrix(x)))
  # the number of samples in the matrix the number of non-genotyped
  # samples remove the loci that are completely missing
  if (length(loci.na) > 0) {
    nind <-
      mean(nrow(as.matrix(x)) -
             colSums(is.na(as.matrix(x)))[-loci.na])
    # the number of samples in the matrix the number of
    # non-genotyped samples
  } else {
    nind <- mean(nrow(as.matrix(x)) - colSums(is.na(as.matrix(x))))
  }
  
  return(nind)
}

# Observed heterozygosity in random subsamples of individuals, by population
# (Schmidt et al. 2021); used by the heterozygosity reports when
# subsample.pop = TRUE.
# x: genlight/dartR object with populations. n.limit: populations with
# fewer individuals are skipped. subsamples: subsample sizes; each
# population uses only the sizes it can supply (at most its number of
# individuals).
# Returns a data.table with one row per population and subsample size:
# res.mean (mean Ho over 10 replicate subsamples), res_SE (standard error
# across those replicates), pop and subsample.
utils.subsample.pop <- function(x,
                                n.limit,
                                subsamples = c(10, 5, 4, 3, 2)){

  x.pops <- seppop(x)
  x.k <- lapply(x.pops, function(y) as.numeric(ploidy(y)))
  x.pops <- lapply(x.pops,as.matrix)

  # Populations below n.limit are skipped, as documented (previously an
  # NA placeholder was stored for them, which data.table::rbindlist
  # rejects, crashing the run for any dataset with a small population)
  keep <- vapply(x.pops, nrow, integer(1)) >= n.limit
  x.pops <- x.pops[keep]
  x.k <- x.k[keep]
  if (length(x.pops) == 0) {
    return(data.table::data.table())
  }

  # Sampling is without replacement, so a population can only supply
  # subsamples up to its own size (reachable when n.limit is below the
  # largest subsample size); populations that can supply none are skipped
  sizes <- lapply(x.pops, function(m) subsamples[subsamples <= nrow(m)])
  keep <- lengths(sizes) > 0
  x.pops <- x.pops[keep]
  x.k <- x.k[keep]
  sizes <- sizes[keep]
  if (length(x.pops) == 0) {
    return(data.table::data.table())
  }

  pops.list <- mapply(function(pop.tmp, k, sz){
    lapply(sz, function(y){
      het_rep(mat = pop.tmp ,samples = y , reps = 10, ploidy = k)
    })
  }, x.pops, x.k, sizes, SIMPLIFY = FALSE)
  pops.list <- lapply(pops.list,data.table::rbindlist)
  pops.list <- lapply(seq_along(pops.list),function(z){
    ptmp <- pops.list[[z]]
    ptmp$pop <- names(pops.list)[z]
    ptmp$subsample <- sizes[[z]]
    return(ptmp)
  })

  return(data.table::rbindlist(pops.list))
}

# Mean and standard error of Ho over `reps` random subsamples of `samples`
# rows (individuals) of the genotype matrix `mat`, drawn without
# replacement. ploidy: one value per row of mat, or a single value. For
# ploidy other than 2, Ho is the gametic heterozygosity d(k - d) /
# choose(k, 2), as in pop.het_fun.
het_rep <- function(mat,samples,reps, ploidy = 2){
  k_all <- rep_len(ploidy, nrow(mat))
  res_tmp <-
    replicate(n = reps, {
      rows <- sample(x = 1:nrow(mat),
                     size = samples,
                     replace = FALSE)
      if (all(ploidy == 2)) {
        het <- mat[rows, , drop = FALSE] == 1
      } else {
        sub <- mat[rows, , drop = FALSE]
        k <- k_all[rows]
        het <- sub * (k - sub) / choose(k, 2)
      }
      mean(colMeans(het, na.rm = TRUE), na.rm = TRUE)
    })
  
  return(data.frame(res.mean = mean(res_tmp) , res_SE = std.error(res_tmp) ))
}
