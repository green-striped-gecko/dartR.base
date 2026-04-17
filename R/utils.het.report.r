# Het calculations
pop.het_fun <- function(df,
                        n.invariant,
                        aHet,
                        bootstrap=TRUE) {
  # rm loci that are all NA
  # otherwise these loci get Ho=0 which would not be correct
  loc.allNA <- colSums(is.na(df)) == nrow(df)
  df <- df[, !loc.allNA]
  
  Ho.loc <- colMeans(df == 1, na.rm = TRUE)
  n_loc.sample <- apply(df, 1, function(y) {
    sum(!is.na(y))
  })
  n_loc <- ncol(df)
  q_freq <- colMeans(df, na.rm = TRUE) / 2
  p_freq <- 1 - q_freq
  He.loc <- 2 * p_freq * q_freq
  n_ind.loc <- apply(df, 2, function(y) {
    sum(!is.na(y))
  })
  ### CP ### Unbiased He (i.e. corrected for sample size) 
  # hard coded for diploid
  uHe.loc <- (2 * as.numeric(n_ind.loc) / (2 * as.numeric(n_ind.loc) - 1)) * He.loc
  
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


# bootstrapping function
pop.het <- function(df,
                    indices,
                    n.invariant = 0,
                    boot_method = "loc",
                    aHet=FALSE) {
  
  df <- df[indices,]
  
  if(boot_method == "loc"){
    df <- t(df)
  }

  res <- pop.het_fun(df,
                     n.invariant = n.invariant,
                     aHet = aHet)
  
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

# standard error function
std.error <- function(x) {
  res <- sd(x, na.rm = TRUE) / sqrt(length(x))
  return(res)
}

utils.subsample.pop <- function(x,
                                n.limit,
                                subsamples = c(10, 5, 4, 3, 2)){
  
  x.pops <- seppop(x)
  x.pops <- lapply(x.pops,as.matrix)
  pops.list <- as.list(1:length(x.pops))
  
for(i in 1:length(x.pops)){
  pop.tmp <- x.pops[[i]]
  if(nrow(pop.tmp) < n.limit){
    pops.list[[i]] <- NA
  }else{
      pops.list[[i]] <- lapply(subsamples, function(y){
        res_tmp <- het_rep(mat = pop.tmp ,samples = y , reps = 10)
      })
  } 
}
  names(pops.list) <- popNames(x)
  pops.list <- lapply(pops.list,data.table::rbindlist)
  pops.list <- lapply(1:nPop(x),function(z){
    ptmp <- pops.list[[z]]
    ptmp$pop <- names(pops.list)[z]
    ptmp$subsample <- subsamples
    return(ptmp)
  })
  
  return(data.table::rbindlist(pops.list))
}

het_rep <- function(mat,samples,reps){
  res_tmp <-
    replicate(n = reps,
              mean(
                colMeans(
                  mat[sample(x = 1:nrow(mat), 
                             size = samples, 
                             replace = FALSE),]== 1,
                  na.rm = TRUE),
                na.rm = TRUE)
    )
  
  return(data.frame(res.mean = mean(res_tmp) , res_SE = std.error(res_tmp) ))
}
