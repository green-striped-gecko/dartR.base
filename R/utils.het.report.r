# bootstrapping function
pop.het <- function(x,
                    indices,
                    n.invariant,
                    aHet=FALSE) {
  
  pop.het_fun <- function(df,
                          n.invariant,
                          aHet) {
    Ho.loc <- colMeans(df == 1, na.rm = TRUE)
    n_loc <- apply(df, 2, function(y) {
      sum(!is.na(y))
    })
    Ho.adj.loc <- Ho.loc * n_loc / (n_loc + n.invariant)
    q_freq <- colMeans(df, na.rm = TRUE) / 2
    p_freq <- 1 - q_freq
    He.loc <- 2 * p_freq * q_freq
    n_ind <- apply(df, 2, function(y) {
      sum(!is.na(y))
    })
    ### CP ### Unbiased He (i.e. corrected for sample size) hard
    # coded for diploid
    uHe.loc <-
      (2 * as.numeric(n_ind) / (2 * as.numeric(n_ind) - 1)) * He.loc
    Hexp.adj.loc <- He.loc * n_loc / (n_loc + n.invariant)
    
    FIS.loc <- 1 - (Ho.loc / He.loc)
    
    if(aHet) {
      all.res <- c(
        Ho.adj.loc = mean(Ho.adj.loc, na.rm = TRUE),
        Hexp.adj.loc = mean(Hexp.adj.loc, na.rm = TRUE)
      )
    } else {
      all.res <- c(
        Ho.loc = mean(Ho.loc, na.rm = TRUE),
        He.loc = mean(He.loc, na.rm = TRUE),
        uHe.loc = mean(uHe.loc, na.rm = TRUE),
        FIS.loc = mean(FIS.loc, na.rm = TRUE)
      )
    }
    
    return(all.res)
  }
  
  df <- x[, indices]
  
  res <- pop.het_fun(df,
                     n.invariant = n.invariant,
                     aHet=aHet)
  
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
