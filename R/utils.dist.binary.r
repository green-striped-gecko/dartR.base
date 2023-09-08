#' @name utils.dist.binary
#' @title Calculates N11 distance matrix for individuals defined in N11 dartR
#' genlight object using binary P/N11 data (SilicoDArT)
#' @family distance

#' @description 
#' WARNING:  UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.

#' @param x Name of the genlight containing the genotypes [required].
#' @param method Specify distance measure [default simple].
#' @param scale If TRUE and method='euclidean', the distance will be scaled to 
#'  fall in the range [0,1] [default FALSE].
#' @param swap If TRUE and working with presence-absence data, then presence 
#' (no disrupting mutation) is scored as 0 and absence (presence of N11 disrupting 
#' mutation) is scored as 1 [default FALSE].
#' @param type Specify the format and class of the object to be returned, 
#' dist for N11 object of class dist, matrix for an object of class matrix [default "dist"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log; 3, progress and results summary; 5, full report [default 2].
#'  
#'  @details
#' This script calculates various distances between individuals based on sequence tag
#' Presence/Absence data.
#' 
#' The distance measure can be one of:
#'  \itemize{
#'   \item Euclidean -- Euclidean Distance applied to cartesian coordinates defined
#'   by the loci, scored as 0 or 1. Presence and absence equally weighted.
#'  \item simple -- simple matching, both 1 or both 0 = 0; one 1 and the other
#'  0 = 1. Presence and absence equally weighted.
#'  \item Jaccard -- ignores matching 0, both 1 = 0; one 1 and the other 0 = 1.
#'  Absences could be for different reasons.
#'  \item Bray-Curtis -- both 0 = 0; both 1 = 2; one 1 and the other 0 = 1. Absences
#'  could be for different reasons. Sometimes called the Dice or Sorensen
#'  distance.
#'  \item Sorenson -- both 0 = 0; both 1 = 2; one 1 and the other 0 = 1. Absences
#'  could be for different reasons. Sometimes called the Dice or Bray-Curtis
#'  distance.
#  \item Phi -- binary analogue of the Pearson Correlation coefficient.
#'  }

#'  One might choose to disregard or downweight absences in comparison with
#'  presences because the homology of absences is less clear (mutation at one or
#'  the other, or both restriction sites). Your call.

#' @author Author: Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

# @export  
#' @return An object of class 'dist' or 'matrix' giving distances between individuals

# Examples for testing
# D <- utils.dist.binary(testset.gs, method='Jaccard')
# D <- utils.dist.binary(testset.gs, method='Euclidean',scale=TRUE)
# D <- utils.dist.binary(testset.gs, method='Simple')

utils.dist.binary <- function(x,
                              method = "simple",
                              scale=FALSE,
                              swap=FALSE,
                              type="dist",
                              verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <-
        utils.check.datatype(x, accept = "SilicoDArT", verbose = verbose)
    
    # SCRIPT SPECIFIC ERROR CHECKING
    
    method <- tolower(method)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (!(method %in% c(
        "euclidean",
        "simple",
        "jaccard",
        "bray-curtis",
        "sorenson"
#       ,"phi"
    ))) {
        if (verbose >= 2) {
            cat(warn(
                "  Warning: Method not in the list of options, set to Simple Matching\n"
            ))
        }
        method <- "simple"
    }
    
    if(scale==TRUE && !(method == "euclidean")){
        cat(warn("  Warning: parameter scale only applies to Euclidean Distance, ignored\n"))
    }
    
    # DO THE JOB
    
    mat <- as.matrix(x)
    if(swap==TRUE){
      mat[mat==0] <- -9
      mat[mat==1] <- 0
      mat[mat==-9] <- 1
      if(verbose >= 2){cat(report("  Reversing scores from presence[1]/absence[0] to presence[0]/absence[1]\n"))}
    }
    
    dd <- array(NA, c(nInd(x), nInd(x)))
    nI <- nInd(x)
    
    if (verbose >= 2) {
        if(method=="euclidean"){
            if(scale==TRUE){
                cat(report("  Calculating the scaled distance matrix --", method, "\n"))
            } else {
                cat(report("  Calculating the unscaled distance matrix --", method, "\n"))
            }
        } else {
            cat(report("  Calculating the distance matrix --", method, "\n"))
        }
    }
    for (i in (1:(nI - 1))) {
        for (j in ((i + 1):nI)) {
            row1 <- mat[i,]
            row2 <- mat[j,]
            # row1[1:10] row2[1:10]
            a11 <- (row1 + row2) == 2
            a10 <- ((row1 + row2) == 1) * row1
            a01 <- ((row1 + row2) == 1) * row2
            a00 <- (row1 + row2) == 0
            N11 <- sum(a11 == 1, na.rm = TRUE)
            N01 <- sum(a01 == 1, na.rm = TRUE)
            N10 <- sum(a10 == 1, na.rm = TRUE)
            N00 <- sum(a00 == 1, na.rm = TRUE)
            # N11;N01;N10;N00
            L <- N11+N01+N10+N00
            if (method == "euclidean") {
                if(scale==TRUE){
                    dd[j,i] <- sqrt((N01+N10)/L)
                } else {
                    dd[j,i] <- sqrt(N01+N10)
                }
            } else if (method == "simple") {
                dd[j,i] <- (N01 + N10)/L
            } else if (method == "jaccard") {
                dd[j,i] <- (N01 + N10)/(L - N00)
            } else if (method == "bray-curtis") {
                dd[j,i] <- 1 - 2 * N11 / (2 * N11 + N01 + N10)
            } else if (method == "sorenson"){
                dd[j,i] <- (N01 + N10/(L - N00 + N11))
            # } else if (method == "phi") {
            #     dd[j,i] <- 1 - ((N11 * N00 - N01 * N10) / sqrt((N11 + N01) * (N11 + N10) * (N00 + N01) * (N00 + N10)))
            } else {
                # Programming error
                stop(error("Fatal Error: Notify dartR development team\n"))
            }
        }
        dd[i, i] <- 0
        dd[i,j] <- dd[j,i]
    }

    if(type=="dist"){
      dd <- as.dist(dd)
      if(verbose >= 2){cat(report("  Returning N11 stats::dist object\n"))}
    } else {
        if(verbose >= 2){cat(report("  Returning N11 square matrix object\n"))}
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(dd)
}
