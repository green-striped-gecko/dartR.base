#' @name gl.join
# Preliminaries -- parameter definitions -----------------
#' @title Combines two dartR genlight objects
#' @family data manipulation
#' 
#' @description
#' This function combines two genlight objects and their associated metadata.
#' The history associated with the two genlight objects is cleared from the new
#' genlight object. Either the individuals/samples must be the same in each genlight
#' object, in which case the new genlight object has the same individuals but combined loci,
#' or the number of loci must be the same in each genlight object in which case the new
#' genlight object has the same loci but combined individuals/samples.

#' The function is typically used to combine datasets from the same service
#' where the files have been split because of size limitations. The data is read
#' in from multiple csv files, then the resultant genlight objects are combined.
#' 
#' This function works with both SNP and Tag P/A data.
#' 
#' @param x1 Name of the first genlight object [required].
#' @param x2 Name of the second genlight object [required].
#' @param method Legacy parameter, issue warning [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#'
#' @details
#' This script joins two genlight objects together along with the associated metadata. if method='sidebyside' (the default), the individuals 
#' in the two genlight objects must be the same and in the same order. The loci are combined.
#' 
#' If method='end2end', the loci in the two genlight objects must be the same and in the same order. The data for the two sets of individuals
#' are combined. Note that if two individuals have the same names, they will be made unique.#' 

#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' # Joining by loci in common, both datasets have the same loci in the same order
#' x1 <- testset.gl[1:7, ]
#' nInd(x1)
#' x2 <- testset.gl[11:14, ]
#' nInd(x2)
#' gl <- gl.join(x1, x2, verbose = 4)
#' nInd(gl)
#' # Joining by individuals in common, both datasets have the same individuals
#' # in the same order
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' x1 <- platypus.gl[, 1:100]
#' nLoc(x1)
#' x2 <- platypus.gl[, 101:200]
#' nLoc(x2)
#' gl <- gl.join(x1, x2, verbose=3)
#' nLoc(gl)
#' 
#' # Join by adding individuals with a set of common loci
#' nInd(testset.gl)
#' x1 <- gl.drop.ind(testset.gl,ind.list=c("AA010915","UC_00126","AA032760","AA013214",
#' "AA011723","AA012411","AA019237","AA019238","AA019239","AA019235","AA019240",
#' "AA019241","AA019242","AA019243"))
#' nInd(x1)
#' x2 <- gl.keep.ind(testset.gl,ind.list=c("AA010915","UC_00126","AA032760","AA013214",
#' "AA011723","AA012411","AA019237","AA019238","AA019239","AA019235","AA019240",
#' "AA019241","AA019242","AA019243"))
#' nInd(x2)
#' gl <- gl.join(x1, x2, verbose=3)
#' nInd(gl)
#' @export
#' @return A new genlight object

gl.join <- function(x1,
                    x2,
                    method=NULL,
                    verbose = NULL) {
  # Preliminaries -------------------
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbose = verbose)
  
  # STANDARD ERROR CHECKING
  
  # CHECK DATATYPE
  datatype1 <- utils.check.datatype(x1, verbose = 0)
  datatype2 <- utils.check.datatype(x2, verbose = 0)
  
  if (!is(x1, "dartR")) {
    class(x1) <- "dartR"
    if (verbose > 2) {
      cat(
        warn(
          "Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"
        )
      )
      cat(
        warn(
          "                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR,
                 please use function dartR2gl\n"
        )
      )
    }
  }
  if (!is(x2, "dartR")) {
    class(x2) <- "dartR"
    if (verbose > 2) {
      cat(
        warn(
          "Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"
        )
      )
      cat(
        warn(
          "                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR,
                 please use function dartR2gl\n"
        )
      )
    }
  }
    if (datatype1 == "SilicoDArT" && datatype2 == "SilicoDArT") {
        if (verbose >= 2) {
            cat(report("  Processing Presence/Absence (SilicoDArT) data\n"))
        }
    } else if (datatype1 == "SNP" && datatype2 == "SNP") {
        if (verbose >= 2) {
            cat(report("  Processing SNP data \n"))
        }
    }
   
    # SCRIPT SPECIFIC ERROR CHECKING
    
    if(!is.null(method)){
      cat(warn("  Warning: The parameter method is deprecated, no longer required"))
      if (method=="join.by.loc"){
        if(verbose >= 3){
          cat(report(" Joining two genlight datasets with the same loci but different individuals\n"))
        }
        flag <- "loc"
      } else if (method=="join.by.ind"){
        if(verbose >= 3){
          cat(report(" Joining two genlight datasets with the same individuals but different loci\n"))
        }
        flag <- "ind"
      } else {
        cat(error("Fatal Error: method parameter is deprecated, no longer required. Please remove from function call\n"))
        stop()
      }
    }
  
  
  if (is.null(method)) {
    
    if (identical(indNames(x1), indNames(x2))) {
      if(verbose >= 3){
        cat(report(" Joining two genlight datasets with the same individuals but different loci\n"))
      }
        flag <- "ind"
    } else if(identical(locNames(x1), locNames(x2) )) {
        if(verbose >= 3){
          cat(report(" Joining two genlight datasets with the same loci but different individuals\n"))
        }
        flag <- "loc"
    } else {
      cat(error("Fatal Error: Individuals or loci in the two files do not match\n"))
      stop()
    }
}
# DO THE JOB --------------
    
if (verbose >= 2) {
  if(flag == "ind"){
    cat(
      report(
        "  Concatenating two genlight objects,",
        substitute(x1),
        "and",
        substitute(x2),
        "with the same individuals, different loci\n"
      )
    )
  }
  
  if(flag == "loc"){
    cat(
      report(
        "  Concatenating two genlight objects,",
        substitute(x1),
        "and",
        substitute(x2),
        "with the same loci, different individuals\n"
      )
    )
  }
}

if(flag == "ind"){
  if (verbose >= 3) {
    cat("    Number of individuals:", nInd(x1), "\n")
    cat("    First genlight object",
        substitute(x1),
        "has",
        nLoc(x1),
        "loci\n")
    cat("    Second genlight object",
        substitute(x2),
        "has",
        nLoc(x2),
        "loci\n")
  }
  
  # Join the two genlight objects
  x <- cbind(x1, x2)
  
  # Join the locus metrics, if they exist
  if (verbose >= 2) {
    cat(report("  Concatenating the locus metrics\n"))
  }
  if (!is.null(x1@other$loc.metrics) &
      !is.null(x2@other$loc.metrics)) {
    x@other$loc.metrics <-
      rbind(x1@other$loc.metrics, x2@other$loc.metrics)
  } else {
    cat(
      warn(
        "  Warning: Input genlight objects lacks locus metrics\n"
      )
    )
  }
  
  # Cater for some locus names being the same in both genlight objects
  if(length(unique(locNames(x))) < length(locNames(x))){
    locNames(x) <- make.unique(locNames(x))
    if(verbose>=3){cat(warn("  Warning: Some locus names in combined genlight object are the same, made unique\n"))}
  }
  
  # Add the ind metrics, assuming they are the same in both genlight objects
  if (verbose >= 2) {
    cat(report("  Adding the individual metrics\n"))
  }
  
  if (!is.null(x1@other$ind.metrics)) {
    x@other$ind.metrics <- x1@other$ind.metrics
  } else if (!is.null(x2@other$ind.metrics)) {
    x@other$ind.metrics <- x2@other$ind.metrics
  } else {
    cat(
      warn(
        "  Warning: Input genlight objects lack individual metrics\n"
      )
    )
  }
 
  # Add the loc metrics flags, set to 1 only if 1 in both genlight objects
  if (verbose >= 2) {
    cat(report("  Setting the locus metrics flags\n"))
  }
  if (!is.null(x1@other$loc.metrics.flags) &
      !is.null(x2@other$loc.metrics.flags)) {
    x@other$loc.metrics.flags$AvgPIC <-
      x1@other$loc.metrics.flags$AvgPIC * x2@other$loc.metrics.flags$AvgPIC
    x@other$loc.metrics.flags$OneRatioRef <-
      x1@other$loc.metrics.flags$OneRatioRef * x2@other$loc.metrics.flags$OneRatioRef
    x@other$loc.metrics.flags$OneRatioSnp <-
      x1@other$loc.metrics.flags$OneRatioSnp * x2@other$loc.metrics.flags$OneRatioSnp
    x@other$loc.metrics.flags$PICRef <-
      x1@other$loc.metrics.flags$PICRef * x2@other$loc.metrics.flags$PICRef
    x@other$loc.metrics.flags$PICSnp <-
      x1@other$loc.metrics.flags$PICSnp * x2@other$loc.metrics.flags$PICSnp
    x@other$loc.metrics.flags$CallRate <-
      x1@other$loc.metrics.flags$CallRate * x2@other$loc.metrics.flags$CallRate
    x@other$loc.metrics.flags$maf <-
      x1@other$loc.metrics.flags$maf * x2@other$loc.metrics.flags$maf
    x@other$loc.metrics.flags$FreqHets <-
      x1@other$loc.metrics.flags$FreqHets * x2@other$loc.metrics.flags$FreqHets
    x@other$loc.metrics.flags$FreqHomRef <-
      x1@other$loc.metrics.flags$FreqHomRef * x2@other$loc.metrics.flags$FreqHomRef
    x@other$loc.metrics.flags$FreqHomSnp <-
      x1@other$loc.metrics.flags$FreqHomSnp * x2@other$loc.metrics.flags$FreqHomSnp
    x@other$loc.metrics.flags$monomorphs <-
      x1@other$loc.metrics.flags$monomorphs * x2@other$loc.metrics.flags$monomorphs
    x@other$loc.metrics.flags$OneRatio <-
      x1@other$loc.metrics.flags$OneRatio * x2@other$loc.metrics.flags$OneRatio
    x@other$loc.metrics.flags$PIC <-
      x1@other$loc.metrics.flags$PIC * x2@other$loc.metrics.flags$PIC
  } else {
    cat(
      warn(
        "  Warning: Input genlight objects and/or output genlight object lacks metrics flags. Flags set to zero\n"
      )
    )
    x@other$loc.metrics.flags$AvgPIC <- 0
    x@other$loc.metrics.flags$OneRatioRef <- 0
    x@other$loc.metrics.flags$OneRatioSnp <- 0
    x@other$loc.metrics.flags$PICRef <- 0
    x@other$loc.metrics.flags$PICSnp <- 0
    x@other$loc.metrics.flags$CallRate <- 0
    x@other$loc.metrics.flags$maf <- 0
    x@other$loc.metrics.flags$FreqHets <- 0
    x@other$loc.metrics.flags$FreqHomRef <- 0
    x@other$loc.metrics.flags$FreqHomSnp <- 0
    x@other$loc.metrics.flags$monomorphs <- 0
    x@other$loc.metrics.flags$OneRatio <- 0
    x@other$loc.metrics.flags$PIC <- 0
  }
}

if(flag=="loc"){
  if (verbose >= 3) {
    cat("    Number of loci:", nLoc(x1), "\n")
    cat("    First genlight object",
        substitute(x1),
        "has",
        nInd(x1),
        "individuals\n")
    cat("    Second genlight object",
        substitute(x2),
        "has",
        nInd(x2),
        "individuals\n")
  }
  
  # Join the two genlight objects
  x <- rbind(x1, x2)
  
  # Join the locus metrics, if they exist
  if (verbose >= 2) {
    cat(report("  Concatenating the individual metrics\n"))
  }
  if (!is.null(x1@other$ind.metrics) &
      !is.null(x2@other$ind.metrics)) {
    x@other$ind.metrics <- dplyr::bind_rows(x1@other$ind.metrics, x2@other$ind.metrics)
  } else {
    cat(
      warn(
        "  Warning: Input genlight objects both lack individual metrics\n"
      )
    )
  }
  
  # Cater for some individual names being the same in both genlight objects
  if(length(unique(indNames(x))) < length(indNames(x))){
    indNames(x) <- make.unique(indNames(x))
    x@other$ind.metrics$id <- indNames(x)
    if(verbose>=3){cat(warn("  Warning: Some individual names in combined genlight object are the same, made unique\n"))}
  }
  
  # Add the locus metrics, assuming they are the same in both genlight objects
  if (verbose >= 2) {
    cat(report("  Adding the locus metrics\n"))
  }
  if (!is.null(x1@other$loc.metrics)) {
    x@other$loc.metrics <- x1@other$loc.metrics
  } else if (!is.null(x2@other$loc.metrics)) {
    x@other$loc.metrics <- x2@other$loc.metrics
  } else {
    cat(
      warn(
        "  Warning: Input genlight objects lack locus metrics\n"
      )
    )
  }
  
  # Add the loc metrics flags, set to 1 only if 1 in both genlight objects
  if (verbose >= 2) {
    cat(report("  Setting the locus metrics flags\n"))
  }
  if (!is.null(x1@other$loc.metrics.flags) &
      !is.null(x2@other$loc.metrics.flags)) {
    x@other$loc.metrics.flags$AvgPIC <-
      x1@other$loc.metrics.flags$AvgPIC * x2@other$loc.metrics.flags$AvgPIC
    x@other$loc.metrics.flags$OneRatioRef <-
      x1@other$loc.metrics.flags$OneRatioRef * x2@other$loc.metrics.flags$OneRatioRef
    x@other$loc.metrics.flags$OneRatioSnp <-
      x1@other$loc.metrics.flags$OneRatioSnp * x2@other$loc.metrics.flags$OneRatioSnp
    x@other$loc.metrics.flags$PICRef <-
      x1@other$loc.metrics.flags$PICRef * x2@other$loc.metrics.flags$PICRef
    x@other$loc.metrics.flags$PICSnp <-
      x1@other$loc.metrics.flags$PICSnp * x2@other$loc.metrics.flags$PICSnp
    x@other$loc.metrics.flags$CallRate <-
      x1@other$loc.metrics.flags$CallRate * x2@other$loc.metrics.flags$CallRate
    x@other$loc.metrics.flags$maf <-
      x1@other$loc.metrics.flags$maf * x2@other$loc.metrics.flags$maf
    x@other$loc.metrics.flags$FreqHets <-
      x1@other$loc.metrics.flags$FreqHets * x2@other$loc.metrics.flags$FreqHets
    x@other$loc.metrics.flags$FreqHomRef <-
      x1@other$loc.metrics.flags$FreqHomRef * x2@other$loc.metrics.flags$FreqHomRef
    x@other$loc.metrics.flags$FreqHomSnp <-
      x1@other$loc.metrics.flags$FreqHomSnp * x2@other$loc.metrics.flags$FreqHomSnp
    x@other$loc.metrics.flags$monomorphs <-
      x1@other$loc.metrics.flags$monomorphs * x2@other$loc.metrics.flags$monomorphs
    x@other$loc.metrics.flags$OneRatio <-
      x1@other$loc.metrics.flags$OneRatio * x2@other$loc.metrics.flags$OneRatio
    x@other$loc.metrics.flags$PIC <-
      x1@other$loc.metrics.flags$PIC * x2@other$loc.metrics.flags$PIC
  } else {
    cat(
      warn(
        "  Warning: Input genlight objectlacks metrics flags. Flags set to zero\n"
      )
    )
    x@other$loc.metrics.flags$AvgPIC <- 0
    x@other$loc.metrics.flags$OneRatioRef <- 0
    x@other$loc.metrics.flags$OneRatioSnp <- 0
    x@other$loc.metrics.flags$PICRef <- 0
    x@other$loc.metrics.flags$PICSnp <- 0
    x@other$loc.metrics.flags$CallRate <- 0
    x@other$loc.metrics.flags$maf <- 0
    x@other$loc.metrics.flags$FreqHets <- 0
    x@other$loc.metrics.flags$FreqHomRef <- 0
    x@other$loc.metrics.flags$FreqHomSnp <- 0
    x@other$loc.metrics.flags$monomorphs <- 0
    x@other$loc.metrics.flags$OneRatio <- 0
    x@other$loc.metrics.flags$PIC <- 0
  }
}

    # Create the history repository, taking the base from X1 if it exists
    if (verbose >= 2) {
      cat(report("  Adding the locus metrics\n"))
    }
    if (!is.null(x1@other$loc.metrics)) {
      x@other$loc.metrics <- x1@other$loc.metrics
    } else if (!is.null(x2@other$loc.metrics)) {
      x@other$loc.metrics <- x2@other$loc.metrics
    } else {
      cat(warn("  Warning: Input genlight objects lack locus metrics\n"))
    }
    
    if (verbose >= 3) {
        cat("    Combined genlight object has", nInd(x), "individuals\n")
        cat("    Combined genlight object has", nLoc(x), "loci\n")
    }
    if (!is.null(x1@other$loc.metrics.flags) &
        !is.null(x2@other$loc.metrics.flags)) {
      x@other$loc.metrics.flags$AvgPIC <-
        x1@other$loc.metrics.flags$AvgPIC * x2@other$loc.metrics.flags$AvgPIC
      x@other$loc.metrics.flags$OneRatioRef <-
        x1@other$loc.metrics.flags$OneRatioRef * x2@other$loc.metrics.flags$OneRatioRef
      x@other$loc.metrics.flags$OneRatioSnp <-
        x1@other$loc.metrics.flags$OneRatioSnp * x2@other$loc.metrics.flags$OneRatioSnp
      x@other$loc.metrics.flags$PICRef <-
        x1@other$loc.metrics.flags$PICRef * x2@other$loc.metrics.flags$PICRef
      x@other$loc.metrics.flags$PICSnp <-
        x1@other$loc.metrics.flags$PICSnp * x2@other$loc.metrics.flags$PICSnp
      x@other$loc.metrics.flags$CallRate <-
        x1@other$loc.metrics.flags$CallRate * x2@other$loc.metrics.flags$CallRate
      x@other$loc.metrics.flags$maf <-
        x1@other$loc.metrics.flags$maf * x2@other$loc.metrics.flags$maf
      x@other$loc.metrics.flags$FreqHets <-
        x1@other$loc.metrics.flags$FreqHets * x2@other$loc.metrics.flags$FreqHets
      x@other$loc.metrics.flags$FreqHomRef <-
        x1@other$loc.metrics.flags$FreqHomRef * x2@other$loc.metrics.flags$FreqHomRef
      x@other$loc.metrics.flags$FreqHomSnp <-
        x1@other$loc.metrics.flags$FreqHomSnp * x2@other$loc.metrics.flags$FreqHomSnp
      x@other$loc.metrics.flags$monomorphs <-
        x1@other$loc.metrics.flags$monomorphs * x2@other$loc.metrics.flags$monomorphs
      x@other$loc.metrics.flags$OneRatioSnp <-
        x1@other$loc.metrics.flags$OneRatioSnp * x2@other$loc.metrics.flags$OneRatioSnp
      x@other$loc.metrics.flags$OneRatioRef <-
        x1@other$loc.metrics.flags$OneRatioRef * x2@other$loc.metrics.flags$OneRatioRef
      x@other$loc.metrics.flags$AvgPIC <-
        x1@other$loc.metrics.flags$AvgPIC * x2@other$loc.metrics.flags$AvgPIC
    } else {
      cat(warn(
        "  Warning: Input genlight objectlacks metrics flags. Flags set to zero\n"
      ))
      x@other$loc.metrics.flags$AvgPIC <- 0
      x@other$loc.metrics.flags$OneRatioRef <- 0
      x@other$loc.metrics.flags$OneRatioSnp <- 0
      x@other$loc.metrics.flags$PICRef <- 0
      x@other$loc.metrics.flags$PICSnp <- 0
      x@other$loc.metrics.flags$CallRate <- 0
      x@other$loc.metrics.flags$maf <- 0
      x@other$loc.metrics.flags$FreqHets <- 0
      x@other$loc.metrics.flags$FreqHomRef <- 0
      x@other$loc.metrics.flags$FreqHomSnp <- 0
      x@other$loc.metrics.flags$monomorphs <- 0
      x@other$loc.metrics.flags$OneRatio <- 0
      x@other$loc.metrics.flags$PIC <- 0
    }
  
  
  # Create the history repository, taking the base from X1 if it exists
  if (verbose >= 2) {
    cat(report("  Adding the history\n"))
  }
  if (is.null(x@other$history)) {
    x@other$history <- list(match.call())
  } else {
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()
  }
  
  if (verbose >= 3) {
    cat("    Number of individuals:", nInd(x1), "\n")
    cat("    Combined genlight object has", nLoc(x), "loci\n")
  }
  
  # FLAG SCRIPT END ---------------
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  # End block -----------------------
  
  return(x)
}
