#' @name gl.impute
#' @title Imputes missing data
#' @family data manipulation

#' @description
#' This function imputes genotypes on a population-by-population basis, where
#' populations can be considered panmictic, or imputes the state for
#' presence-absence data.
#' 
#' @param x Name of the genlight object containing the SNP or presence-absence
#' data [required].
#' @param method Imputation method, either "frequency" or "HW" or "neighbour" 
#' or "random" or "beagle" [default "neighbour"].
#' @param fill.residual Should any residual missing values remaining after 
#' imputation be set to 0, 1, 2 at random, taking into account global allele 
#' frequencies at the particular locus [default TRUE].
#' @param parallel A logical indicating whether multiple cores -if available-
#' should be used for the computations (TRUE), or not (FALSE); requires the
#' package parallel to be installed [default FALSE].
#' @param beagle.bin.path Path of beagle.27Feb25.75f.jar file [default getwd())].
#' @param plink.bin.path Path of PLINK binary file [default getwd())].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' 
#' @details
#' We recommend that imputation be performed on sampling locations, before
#' any aggregation. Imputation is achieved by replacing missing values using
#' either of five methods:
#' \itemize{
#' \item If "frequency", genotypes scored as missing at a locus in an individual
#'  are imputed using the average allele frequencies at that locus in the 
#'  population from which the individual was drawn.
#' \item If "HW", genotypes scored as missing at a locus in an individual are 
#' imputed by sampling at random assuming Hardy-Weinberg equilibrium. Applies 
#' only to genotype data.
#' \item If "neighbour", substitute the missing values for the focal individual
#'  with the values taken from the nearest neighbour. Repeat with next nearest
#'  and so on until all missing values are replaced.
#' \item if "random", missing data are substituted by random values (0, 1 or 2). 
#' \item if "beagle", missing data is imputed using using BEAGLE 
#' (beagle.27Feb25.75f.jar), which infers missing genotypes by modelling 
#' shared haplotype patterns among individuals. Beagle can be downloaded using 
#' the following link: 
#' 
#' https://faculty.washington.edu/browning/beagle/beagle.html#download.
#' 
#' After downloading the Beagle binary move it to your working directory. 
#' 
#' For method = "beagle" it is also required to download the binary file of PLINK 1.9
#' and move it to your working directory. The binary file can be downloaded from:
#' 
#' \url{https://www.cog-genomics.org/plink/}
#' }
#'   The nearest neighbour is the one at the smallest Euclidean distance from 
#'   the focal individual.
#'   
#'   The advantage of this approach is that it works regardless of how many
#'   individuals are in the population to which the focal individual belongs,
#'   and the displacement of the individual is haphazard as opposed to:
#'   (a) Drawing the individual toward the population centroid (HW and Frequency).
#'   (b) Drawing the individual toward the global centroid (glPCA).
#'   
#' Note that loci that are missing for all individuals in a population are not 
#' imputed with method 'frequency' or 'HW' and can give unpredictable results
#' for particular individuals using 'neighbour'.
#' 
#'  Consider using the function 
#' \code{\link{gl.filter.allna}} with by.pop=TRUE to remove them first.
#'@references
#'\itemize{
#'\item Browning, B. L., & Browning, S. R. (2016). Genotype imputation with 
#'millions of reference samples. American Journal of Human Genetics, 98(1), 
#'116–126. https://doi.org/10.1016/j.ajhg.2015.11.020
#' }
#' @author Custodian: Luis Mijangos 
#' (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'  \donttest{
#' require("dartR.data")
#' # SNP genotype data
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' gl <- gl.filter.callrate(platypus.gl,threshold=0.95)
#' gl <- gl.filter.allna(gl)
#' gl <- gl.impute(gl, method="frequency") 
#' # Sequence Tag presence-absence data
#' gs <- gl.filter.callrate(testset.gs,threshold=0.95)
#' gl <- gl.filter.allna(gl)
#' #gs <- gl.impute(gs, method="neighbour")
#' }
#' gl <- gl.impute(platypus.gl,method ="random")
#' 
#' @export
#' @return A genlight object with the missing data imputed.

gl.impute <-  function(x,
                       method = "neighbour",
                       fill.residual = TRUE,
                       parallel = FALSE,
                       beagle.bin.path = getwd(),
                       plink.bin.path = getwd(),
                       verbose = NULL) {
  
  x_hold <- x
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
    if(method=="neighbor"){
      method <- "neighbour"
    }
  
  # DO THE JOB
  #check type
  
  fbm <- .fbm_or_null(x)
  #separating populations
  
  if (method == "frequency" | method == "HW") {
    pop_list_temp <- seppop(x)
    pop_list <- list()
    
    for (y in pop_list_temp) {
      loci_all_nas <- sum(glNA(y) > nInd(y))
      nas_number <- sum(glNA(y)) / 2
      number_imputations <- nas_number - (loci_all_nas * nInd(y))
      
      if (method == "frequency") {
        if (verbose >= 2){
          cat(report("  Imputation based on average allele frequencies, population-wise\n"))
        }
        if (verbose >= 2 & loci_all_nas >= 1) {
          cat(
            warn(
              "  Warning: Population ",
              popNames(y),
              " has ",
              loci_all_nas,
              " loci with all missing values.\n"
            )
          )
          if (verbose >= 3) {
            cat(
              report(
                "  Method= 'frequency':",
                number_imputations,
                "values to be imputed.\n"
              )
            )
          }
        }
        
        q_allele <- glMean(y)
        pop_matrix <- as.matrix(y)
        loc_na <- which(is.na(pop_matrix), arr.ind = TRUE)
pop_matrix[loc_na] <- unname(unlist(lapply(q_allele[loc_na[, 2]], function(x) {
            return(as.numeric(s_alleles(q_freq = x)))
          })))
        
        if (is.null(fbm)) y@gen <- matrix2gen(pop_matrix, parallel = parallel) else y@fbm[] <-pop_matrix
        pop_list <- c(pop_list, y)
      }
      
      if (method == "HW") {
        if (verbose >= 2){
          cat(report("  Imputation based on average allele HW sampling, population-wise\n"))
        }
        if (verbose >= 2 & loci_all_nas >= 1) {
          cat(
            warn(
              "  Warning: Population ",
              popNames(y),
              " has ",
              loci_all_nas,
              " loci with all missing values.\n"
            )
          )
          if (verbose >= 3) {
            cat(report(
              "  Method= 'HW':",
              number_imputations,
              "values to be imputed.\n"
            ))
          }
        }
        
        q_allele <- glMean(y)
        pop_matrix <- as.matrix(y)
        loc_na <- which(is.na(pop_matrix), arr.ind = TRUE)
        pop_matrix[loc_na] <-
          unname(unlist(lapply(q_allele[loc_na[, 2]], function(x) {
            return(sample_genotype(q_freq = x))
          })))
        
        if(is.null(fbm)) y@gen <- matrix2gen(pop_matrix, parallel = parallel) else y@fbm[] <-pop_matrix
        pop_list <- c(pop_list, y)
      }
    }
    
    # if more than 1 population
    if (length(pop_list) > 1) {
      x3 <- NULL
      # merge back populations
      for (pop in pop_list) {
        x3 <- rbind(x3, pop)
      }
    }
    
    # if 1 population
    if (length(pop_list) == 1) {
      x3 <- pop_list[[1]]
    }
    
  }
  
  if (method == "neighbour") {
    
    if (verbose >= 2) {
      cat(report("  Imputation based on drawing from the nearest neighbour\n"))
    }
    
    ## ---- Optional: per-population diagnostics (does not affect imputation) ----
    if (verbose >= 2) {
      pop_list_temp <- seppop(x)
      
      for (k in seq_along(pop_list_temp)) {
        yy <- pop_list_temp[[k]]
        
        loci_all_nas <- sum(glNA(yy) >= nInd(yy))
        nas_number <- sum(glNA(yy)) / 2
        number_imputations <- nas_number - (loci_all_nas * nInd(yy))
        
        if (loci_all_nas >= 1) {
          pop_name <- names(pop_list_temp)[k]
          if (is.null(pop_name) || pop_name == "") {
            # fallback if list not named
            pop_name <- tryCatch(as.character(unique(pop(yy))[1]), error = function(e) "UNKNOWN")
          }
          
          cat(warn(
            "  Warning: Population ", pop_name,
            " has ", loci_all_nas, " loci with all missing values.\n",
            sep = ""
          ))
          
          if (verbose >= 3) {
            cat(report(
              "  Method = 'neighbour': ", number_imputations,
              " values to be imputed (excluding all-NA loci).\n",
              sep = ""
            ))
          }
        }
      }
    }
    
    ## ---- Main imputation ----
    x3 <- x
    x_matrix <- as.matrix(x)  # numeric matrix, nInd x nLoc, with NA
    
    D <- gl.dist.ind(
      x,
      method = "Euclidean",
      verbose = 0,
      plot.display = FALSE,
      type = "matrix"
    )
    
    # Robust distance handling:
    D <- as.matrix(D)
    D[is.na(D)] <- Inf     # if any undefined distances, make them last
    diag(D) <- Inf         # never pick self as neighbour
    
    n <- nInd(x)
    
    for (i in seq_len(n)) {
      
      miss <- is.na(x_matrix[i, ])
      if (!any(miss)) next
      
      ord <- order(D[i, ], decreasing = FALSE)  # neighbour indices by increasing distance
      
      # Walk neighbours until all missing loci filled (or we run out of neighbours)
      for (j in ord) {
        if (!any(miss)) break
        
        cand <- x_matrix[j, miss]       # values at missing loci from neighbour j
        ok <- !is.na(cand)              # loci neighbour can actually impute
        
        if (any(ok)) {
          miss_pos <- which(miss)
          fill_pos <- miss_pos[ok]
          x_matrix[i, fill_pos] <- cand[ok]
          miss[fill_pos] <- FALSE
        }
      }
      
      if (any(miss) && verbose >= 1) {
        cat(important(
          "  Unable to fully impute individual ", indNames(x)[i],
          " (", sum(miss), " loci remain NA)\n",
          sep = ""
        ))
      }
    }
    
    if (is.null(fbm)) {
      x3@gen <- matrix2gen(x_matrix, parallel = parallel)
    } else {
      x3@fbm[] <- x_matrix
    }
  }
  
  if (method == "random") {
    pop_list_temp <- seppop(x)
    pop_list <- list()
    
    for (y in pop_list_temp) {
      loci_all_nas <- sum(glNA(y) > nInd(y))
      nas_number <- sum(glNA(y)) / 2
      number_imputations <- nas_number - (loci_all_nas * nInd(y))
    }
    
    if (verbose >= 2 & loci_all_nas >= 1) {
      cat(
        warn(
          "  Warning: Population ",
          popNames(y),
          " has ",
          loci_all_nas,
          " loci with all missing values.\n"
        )
      )
      if (verbose >= 3) {
        cat(report(
          "  Method= 'random':",
          number_imputations,
          "values to be imputed.\n"
        ))
      }
    }

    x3 <- x
    
    x_matrix <- as.matrix(x)
    loc_na <- which(is.na(x_matrix), arr.ind = TRUE)
    x_matrix[loc_na] <- sample(c(0:2),size=nrow(loc_na),replace = TRUE)
    
    if (is.null(fbm)) x3@gen <- matrix2gen(x_matrix, parallel = parallel) else x3@fbm[] <-x_matrix
    
  }
  
  if (method == "beagle") {
    # FUNCTION SPECIFIC ERROR CHECKING check if packages is installed
    pkg <- "R.utils"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    # Rename scaffolds with only one SNP, which breaks beagle
    if(length(x$chromosome)>0){
    chr <- as.character(x@chromosome)
    singletons <- names(which(table(chr) == 1L))
    chr[chr %in% singletons] <- ""
    x@chromosome <- factor(chr)
    }
    
    pop_list_temp <- seppop(x)
    pop_list <- list()
    
    for (y in pop_list_temp) {
      loci_all_nas <- sum(glNA(y) > nInd(y))
      nas_number <- sum(glNA(y)) / 2
      number_imputations <- nas_number - (loci_all_nas * nInd(y))
    }
    
    if (verbose >= 2 & loci_all_nas >= 1) {
      cat(
        warn(
          "  Warning: Population ",
          popNames(y),
          " has ",
          loci_all_nas,
          " loci with all missing values.\n"
        )
      )
      if (verbose >= 3) {
        cat(report(
          "  Method= 'beagle':",
          number_imputations,
          "values to be imputed.\n"
        ))
      }
    }
    
    x_tmp <- x
    x_tmp@position <- 1:nLoc(x_tmp)
    gl2vcf(x_tmp,
           plink.bin.path = plink.bin.path,
           outpath = tempdir())
    system(paste0(
      "java -Xmx20g -jar ",beagle.bin.path,"/beagle.27Feb25.75f.jar gt=",
      tempdir(),
      "/gl_vcf.vcf out=",
      tempdir(),
      "/imputed"))
    R.utils::gunzip(paste0(tempdir(),
                  "/imputed.vcf.gz"),
                  overwrite =T)

    
    if (!is.null(fbm)) {
      x3 <- gl.read.vcf(paste0(tempdir(),
                               "/imputed.vcf"),
                        fbm = TRUE,
                        verbose = 0)
    }else{
      x3 <- gl.read.vcf(paste0(tempdir(),
                               "/imputed.vcf"),
                        verbose = 0)
      
    }
    
  }
  
  if(fill.residual==TRUE){
    
    q_allele <- glMean(x3)
    pop_matrix <- as.matrix(x3)
    loc_na <- which(is.na(pop_matrix), arr.ind = TRUE)
    pop_matrix[loc_na] <- unname(unlist(lapply(q_allele[loc_na[, 2]], function(x) {
      return(as.numeric(s_alleles(q_freq = x)))
    })))
    if (is.null(fbm)) x3@gen <- matrix2gen(pop_matrix, parallel = parallel) else x3@fbm[] <-pop_matrix

    if(verbose>=2){
    cat(report("  Residual missing values were filled randomly drawing from the global allele profiles by locus\n"))
    }
  }
  
  x3$chromosome <- x@chromosome
  x3$position <- x$position
  x3$ploidy <- x$ploidy
  x3$strata <- x$strata
  x3$hierarchy <- x$hierarchy
  x3$other <- x$other
  x3$loc.all <- x$loc.all
  pop(x3) <- pop(x)

  x3 <- gl.compliance.check(x3, verbose = 0)
  
  if(verbose>=3){
    
    pop_list_before <- seppop(x_hold)
    all_nas_before <- sum(unlist(lapply(pop_list_before,function(y){
       sum(glNA(y) > nInd(y))
    })))
    x_matrix_before <- as.matrix(x_hold)
    nas_before <- sum(is.na(x_matrix_before))
    
    pop_list_after <- seppop(x3)
    all_nas_after <- sum(unlist(lapply(pop_list_after,function(y){
      sum(glNA(y) > nInd(y))
    })))

    x_matrix_after <- as.matrix(x3)
    nas_after <- sum(is.na(x_matrix_after))
    imputed <- nas_before - nas_after
    
    cat("  Imputation method:",method,"\n")
    cat("  No. of missing values before imputation:",nas_before,"\n")
    cat("  No. of loci with all NA's for any one population before imputation:",all_nas_before,"\n")
    cat("  No. of values imputed:",imputed,"\n")
    cat("  No. of missing values after imputation:",nas_after,"\n")
    cat("  No. of loci with all NA's for any one population after imputation:",all_nas_after,"\n")
  }

  # ADD TO HISTORY
  x3@other$history <- x@other$history
  nh <- length(x3@other$history)
  x3@other$history[[nh + 1]] <- match.call()
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  return(x3)
}
