#' @name gl.cross
#' @title Generates random crosses between fathers and mothers
#' @family simulation
#' 
#' @description
#' Generates random crosses between fathers (in one genlight object) and
#' mothers (in a second genlight object) then randomly selects a specified
#' number of offspring to retain.
#' 
#' @param fathers Genlight object of potential fathers [required].
#' @param mothers Genlight object of potential mothers simulated [required].
#' @param broodsize Number of offspring per mother [required].
#' @param sexratio Sex ratio of simulated offspring [default 0.5].
#' @param n Number of offspring to retain [default 1000 or mothers*broodsize whichever is the lesser]
#' @param error.check If TRUE, will perform error checks on the provided parameters [default TRUE]
#' @param compliance.check If TRUE, will perform a compliance check on the resultant
#' genlight object before returning it [default TRUE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#' 
#' @details
#' This script is to be used in conjunction with gl.subsample.ind() applied
#' initially to a base genlight object containing initial male and female 
#' genotypes. The workflow is to
#' 
#' (a) Select the males from the base genlight object using gl.keep.pop() with 
#' pop.list="male" and the as.pop parameter set to sex.
#' (b) Select the females from the base genlight object using gl.keep.pop() with 
#' pop.list="female" and the as.pop parameter set to sex.
#' (c) Subsample a cohort of males for breeding and a cohort of females for breeding
#' using gl.subsample.ind() and the replace parameter as follows:
#' 
#' To enforce monogamy -- generate the fathers and mothers from the base genlight object
#' using gl.subsample.ind() with replace=FALSE.
#' To admit polygyny -- generate the fathers from the base genlight object using 
#' gl.subsample.ind() with replace=FALSE and the mothers from the base genlight
#' object using gl.subsample.ind() with replace=TRUE.
#' To admit polyandry -- generate the fathers from the base genlight object using 
#' gl.subsample.ind() with replace=TRUE and the mothers from the base genlight
#' object using gl.subsample.ind() with replace=FALSE.
#' To admit promiscuity -- generate the fathers and mothers from the base genlight object using
#' gl.subsample.ind() with replace=TRUE.
#' 
#' These are simple scenarios that leave the number of maternal mates per 
#' father (polygyny) and the number of paternal mates per mother (polyandry) 
#' to chance, depending on the random selection of males and females with 
#' replacement from the base genlight object.
#' 
#' (d) Cross the males with the females using gl.sim.crosses() retaining a
#' subset of offspring at random.
#' 
#' So the input for this function is a genlight object with a sample of male
#' individuals (fathers) selected from a larger set at random with or without 
#' replacement; a similar sample of female individuals in a second genlight 
#' object (mothers); specified broodsize; and desired offspring sex ratio.
#' 
#' Set check.error to FALSE if using this script in simulations
#' 
#' @author Custodian: Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#' #Simulate 10 potential fathers
#' gl.fathers <- glSim(10, 20, ploidy=2)
#' #Simulate 10 potential mothers
#' gl.mothers <- glSim(10, 20, ploidy=2)
#' gl.cross(gl.fathers, gl.mothers, 2, sexratio=0.5)
#' 
#' @importFrom stats runif
#' @export
#' @return A genlight object with n offspring of both sexes.

gl.cross <- function(fathers,
                             mothers,
                             broodsize=10,
                             sexratio = 0.5,
                             n = 1000,
                             error.check = TRUE,
                             compliance.check = TRUE,
                             verbose=NULL) {
  
  
  if(error.check){
    
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)
    
    # CHECK DATATYPE
    cat(report("father --"))
    datatype.dad <- utils.check.datatype(fathers, verbose = verbose)
    cat(report("mother --"))
    datatype.mum <- utils.check.datatype(mothers, verbose = verbose)
    
    # CHECK BROODSIZE
    if(broodsize < 0){
      cat(warn("  Error: Brood Size must be a positive interger\n    Set to 10\n"))
    }
    
    # CHECK SEXRATIO
    if(sexratio < 0 | sexratio > 1){
      cat(warn("  Error: Sex ratio must be in the range 0 to 1\n    Set to 0.5\n"))
    }
    
    # CHECK TOTAL NUMBER OF OFFSPRING
    noff <- nInd(mothers) * broodsize
    if(noff < n){
      cat(warn("  Error: Sum of broods less than specified number of offspring to return, returning",noff,"\n"))
    }
  }
  
  # Generate maternal haplotypes (ova)  
  mmat <- as.matrix(mothers)
  mhet <- sum(mmat == 1,na.rm=TRUE)
  ova <- array(data=NA,dim=dim(mmat))
  hold <- NULL
  for (i in 1:broodsize){
    ova <- ifelse(mmat == 1, sample(c(0, 2), mhet, replace = T), mmat) #ovum from each mother
    ova <- rbind(hold,ova)
    hold <- ova
  }
  
  # Generate paternal haplotypes (sperm)  
  fmat <- as.matrix(fathers)
  fhet <- sum(fmat == 1,na.rm=TRUE)
  sperm <- array(data=NA,dim=dim(fmat))
  hold <- NULL
  for (i in 1:broodsize){
    sperm <- ifelse(fmat == 1, sample(c(0, 2), fhet, replace = T), fmat) #sperm from each father
    sperm <- rbind(hold,sperm)
    hold <- sperm
  }

  # Generate offspring (zygote) genotypes, taking advantage of the dartR coding
  offmat <- (ova + sperm)/2
  
  gl2 <-
    new(
      "genlight",
      gen = offmat,
      ind.names = paste0("Po_", 1:noff),
      loc.names = locNames(mothers),
      ploidy = rep(2, nrow(offmat))
    )

  # set sex ratio
  sr <- factor(ifelse(runif(nInd(gl2)) < sexratio, "female", "male"))
  gl2@other$ind.metrics$sex <- sr
  
  if(compliance.check){
    gl2 <- gl.compliance.check(gl2, verbose=0)
  }  
  
  if(error.check){
    # ADD TO HISTORY
    nh <- length(gl2@other$history)
    gl2@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END ---------------
    
    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }
  }
  
  return(gl2)
}
