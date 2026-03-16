#' @name utils.dist.ind.snp
#' @title Calculates a distance matrix for individuals defined in a 
#' genlight object using SNP data (DArTseq)
#' @family distance
#' 
#' @description 
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED 
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE 
#' OUTCOMES.
#' 
#' @param x Name of the genlight containing the genotypes [required].
#' @param method Specify distance measure [default Euclidean].
#' @param scale If TRUE and method='Euclidean', the distance will be scaled to 
#'  fall in the range [0,1] [default FALSE].
#' @param type Specify the format and class of the object to be returned, 
#' dist for a object of class dist, matrix for an object of class matrix
#'  [default "dist"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log; 3, progress and results summary; 5, full report [default 2].
#'  
#' @details
#' This script calculates various distances between individuals based on 
#' SNP genotypes.
#'
#' The distance measure can be one of:
#'  \itemize{
#'   \item Euclidean -- Euclidean Distance applied to Cartesian coordinates
#'    defined by the loci, scored as 0, 1 or 2. 
#'  \item Simple -- simple mismatch, 0 where no alleles are shared, 1 where one
#'  allele is shared, 2 where both alleles are shared. 
#'  \item Absolute -- absolute mismatch, 0 where no alleles are shared, 1 where
#'  one or both alleles are shared.
#'  \item Czekanowski (or Manhattan) calculates the city block metric distance
#'  by summing the scores on each axis (locus).
#'  }
#'  
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
# @export
#' @return An object of class 'dist' or 'matrix' giving distances between
#'  individuals

# Examples for testing
# D <- utils.dist.ind.snp(testset.gl, method='Manhattan')
# D <- utils.dist.ind.snp(testset.gl, method='Euclidean', scale = TRUE)
# D <- utils.dist.ind.snp(testset.gl, method='Simple')

utils.dist.ind.snp <- function(x,
                               method  = "Euclidean",
                               scale   = FALSE,
                               type    = "dist",
                               verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)
  
  # FUNCTION-SPECIFIC SETTINGS
  method <- tolower(method)
  
  if (!(method %in% c("euclidean", "simple", "absolute", "czekanowski", "manhattan"))) {
    if (verbose >= 2) {
      cat(warn("  Warning: Method not in the list of options, set to Euclidean\n"))
    }
    method <- "euclidean"
  }
  
  if (scale && method != "euclidean") {
    cat(warn("  Warning: parameter scale only applies to Euclidean Distance, ignored\n"))
  }
  
  # DO THE JOB
  mat <- as.matrix(x)
  
  nI <- nInd(x)
  nL <- nLoc(x)
  dd <- array(NA_real_, c(nI, nI))
  
  # if (verbose >= 2) {
    if (method == "euclidean") {
      if (scale) {
        cat(report("  Calculating the scaled distance matrix --", method, "\n"))
      } else {
        cat(report("  Calculating the unscaled distance matrix --", method, "\n"))
      }
    } else {
      cat(report("  Calculating the distance matrix --", method, "\n"))
    }
    
    # DO THE JOB
    
     mat <- as.matrix(x)
    
    # dd <- array(NA, c(nInd(x), nInd(x)))
    # nI <- nInd(x)
    # nL <- nLoc(x)
     
    if (verbose >= 2) {
        if(method=="euclidean"){
            if(scale==TRUE){
                cat(report("  Calculating the scaled distance matrix --", method, "\n"))
            } else {
                cat(report("  Calculating the unscaled distance matrix --", method, "\n"))
            }
        } #else {
          #dd[j, i] <- sqrt(sum(sq, na.rm = TRUE))
        #}
    }
     
    # for (i in (1:(nI - 1))) {
    #     for (j in ((i + 1):nI)) {
    #         row1 <- mat[i,]
    #         row2 <- mat[j,]
    #         
    #         if (method == "euclidean") {
    #           sq <- (row1-row2)**2
    #           sq <- sq[!is.na(sq)]
    #           L <- length(sq)
    #             if(scale==TRUE){
    #                 dd[j,i] <- sqrt(sum(sq)/L)
    #             } else {
    #                 dd[j,i] <- sqrt(sum(sq))
    #             }
    #         } else if (method == "simple") {
    #           row <- array(1,dim=nL)
    #           row[((row1 + row2) == 4)] <- 2
    #           row[((row1 + row2) == 0)] <- 0
    #           row[is.na(row1 + row2)] <- NA
    #           row <- row[!is.na(row)]
    #           L <- length(row)
    #           dd[j,i] <- 1 - sum(row)/(2*L)
    #         } else if (method == "absolute") {
    #           row <- array(1,dim=nL)
    #           row[((row1 + row2) == 0)] <- 0
    #           row[is.na(row1 + row2)] <- NA
    #           row <- row[!is.na(row)]
    #           L <- length(row)
    #           dd[j,i] <- 1 - sum(row)/(L)
    #         } else if (method == "manhattan" || method == "czekanowski") {
    #           modq <- abs(row1-row2)
    #           modq <- modq[!is.na(modq)]
    #           L <- length(modq)
    #           dd[j,i] <- sum(modq)/(2*L)
    #         } else {
    #             # Programming error
    #             stop(error("Fatal Error: Notify dartR development team\n"))
    #         }
    #     }
    #     dd[i, i] <- 0
    #     dd[i,j] <- dd[j,i]
    # }
    
    dist_mod <- function(){}  #to hack package checking...
    Rcpp::cppFunction(
'   
    // Compute a symmetric distance matrix for rows of x
    // method: "euclidean", "simple", "absolute", "manhattan" or "czekanowski"
    // scale: if true, normalize Euclidean distances by number of loci
    NumericMatrix dist_mod(const NumericMatrix& x,
                           std::string method = "euclidean",
                           bool scale = true) {
      // Number of individuals (rows) and number of loci (columns)
      int nI = x.nrow();
      int nL = x.ncol();
      
      // Initialize output distance matrix with zeros
      NumericMatrix dd(nI, nI);
      
      // Loop over all pairs (i,j), only computing j>i, then mirror
      for (int i = 0; i < nI; ++i) {
        dd(i, i) = 0.0;  // distance from self is zero
        for (int j = i + 1; j < nI; ++j) {
          double val = 0.0;
          
          // Euclidean distance: sqrt(sum((x_i - x_j)^2) / L)
          if (method == "euclidean") {
            double sumsq = 0.0;
            int L = 0;  // count of non-missing values
            for (int k = 0; k < nL; ++k) {
              double v1 = x(i, k), v2 = x(j, k);
              if (NumericMatrix::is_na(v1) || NumericMatrix::is_na(v2)) continue;
              double d = v1 - v2;
              sumsq += d * d;
              ++L;
            }
            if (L > 0) {
              // optionally scale by L
              val = scale ? std::sqrt(sumsq / L) : std::sqrt(sumsq);
            } else {
              val = NA_REAL;  // all values missing
            }
            
            // "simple" distance: based on allele sharing rules
          } else if (method == "simple") {
            double sumrow = 0.0;
            int L = 0;
            for (int k = 0; k < nL; ++k) {
              double v1 = x(i, k), v2 = x(j, k);
              if (NumericMatrix::is_na(v1) || NumericMatrix::is_na(v2)) continue;
              double s = v1 + v2;
              double rv = 1.0;
              if (s == 4.0) rv = 2.0;  // both homozygous same
              else if (s == 0.0) rv = 0.0;  // both homozygous opposite
              sumrow += rv;
              ++L;
            }
            // distance = 1 - (sum of scores) / (2 * L)
            val = (L > 0 ? 1.0 - sumrow / (2.0 * L) : NA_REAL);
            
            // "absolute" distance: 1 - proportion of shared alleles
          } else if (method == "absolute") {
            double sumrow = 0.0;
            int L = 0;
            for (int k = 0; k < nL; ++k) {
              double v1 = x(i, k), v2 = x(j, k);
              if (NumericMatrix::is_na(v1) || NumericMatrix::is_na(v2)) continue;
              double s = v1 + v2;
              // if both alleles zero, rv = 0 else rv = 1
              double rv = (s == 0.0 ? 0.0 : 1.0);
              sumrow += rv;
              ++L;
            }
            // distance = 1 - (sumrow / L)
            val = (L > 0 ? 1.0 - sumrow / L : NA_REAL);
            
            // "manhattan" or "czekanowski": normalized L1 distance
          } else if (method == "manhattan" || method == "czekanowski") {
            double sumabs = 0.0;
            int L = 0;
            for (int k = 0; k < nL; ++k) {
              double v1 = x(i, k), v2 = x(j, k);
              if (NumericMatrix::is_na(v1) || NumericMatrix::is_na(v2)) continue;
              sumabs += std::abs(v1 - v2);
              ++L;
            }
            // distance = sumabs / (2 * L)
            val = (L > 0 ? sumabs / (2.0 * L) : NA_REAL);
            
          } else {
            // Unknown method: set as NA
            val = NA_REAL;
          }
          
          // Fill symmetric entries
          dd(i, j) = val;
          dd(j, i) = val;
        }
      }
      
      return dd;
    }'
    )
    
    dd <- dist_mod(mat, method = method, scale = scale)

    if(type=="dist"){
      dd <- as.dist(dd)
      if(verbose >= 2){cat(report("  Returning a stats::dist object\n"))}
    } else {
        if(verbose >= 2){cat(report("  Returning a square matrix object\n"))}
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
  # }
  
  # Fill diagonal and mirror upper triangle
  dd <- as.matrix(dd)
  diag(dd) <- 0
  dd[upper.tri(dd)] <- t(dd)[upper.tri(dd)]
  
  if (type == "dist") {
    dd <- as.dist(dd)
    if (verbose >= 2) cat(report("  Returning a stats::dist object\n"))
  } else {
    if (verbose >= 2) cat(report("  Returning a square matrix object\n"))
  }
  
  # FLAG SCRIPT END
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(dd)
}
