#' @name gl.join
# Preliminaries -- parameter definitions -----------------
#' @title Combines two dartR genlight objects
#' @family data manipulation
#' 
#' @description
#' This function combines two genlight objects and their associated metadata.
#' The history of the new genlight object is taken from the first genlight
#' object, with the gl.join call appended. Either the individuals/samples must
#' be the same in each genlight
#' object, in which case the new genlight object has the same individuals but combined loci,
#' or the loci must be the same in each genlight object in which case the new
#' genlight object has the same loci but combined individuals/samples.

#' The function is typically used to combine datasets from the same service
#' where the files have been split because of size limitations. The data is read
#' in from multiple csv files, then the resultant genlight objects are combined.
#' 
#' This function works with both SNP and Tag P/A data.
#' 
#' @param x1 Name of the first genlight object [required].
#' @param x2 Name of the second genlight object [required].
#' @param method Legacy parameter, deprecated and no longer required.
#' Accepted legacy values for backward compatibility: 'join.by.loc' or
#' 'end2end' (same loci, individuals combined); 'join.by.ind' or
#' 'sidebyside' (same individuals, loci combined) [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#'
#' @details
#' This script joins two genlight objects together along with the associated
#' metadata. The join mode is detected automatically: if the two objects have
#' identical individuals (in the same order), the loci are combined; if they
#' have identical loci (in the same order), the individuals are combined.
#' Individual names duplicated across the two objects are made unique.
#' The legacy method parameter is deprecated and no longer required; the
#' historical values 'join.by.loc'/'end2end' (join by shared loci) and
#' 'join.by.ind'/'sidebyside' (join by shared individuals) remain accepted
#' for backward compatibility, and the requested join is validated against
#' the data.

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
#' @return A new genlight object with the combined data and metadata.
#' @export

gl.join <- function(x1,
                    x2,
                    method=NULL,
                    verbose = NULL) {
  # Preliminaries -------------------
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # Capture the argument expressions for messages (before evaluation)
  x1.name <- paste(deparse(substitute(x1)), collapse = "")
  x2.name <- paste(deparse(substitute(x2)), collapse = "")
  
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
    if (datatype1 != datatype2) {
        stop(error(
            "Fatal Error: The two genlight objects hold different datatypes (",
            datatype1, " and ", datatype2, ") and cannot be joined\n"
        ))
    }
    if (datatype1 == "SilicoDArT") {
        if (verbose >= 2) {
            cat(report("  Processing Presence/Absence (SilicoDArT) data\n"))
        }
    } else if (datatype1 == "SNP") {
        if (verbose >= 2) {
            cat(report("  Processing SNP data \n"))
        }
    }

    # SCRIPT SPECIFIC ERROR CHECKING

    if(!is.null(method)){
      if (verbose >= 2) {
        cat(warn("  Warning: The parameter method is deprecated, no longer required\n"))
      }
      if (method %in% c("join.by.loc", "end2end")){
        if(verbose >= 3){
          cat(report(" Joining two genlight datasets with the same loci but different individuals\n"))
        }
        flag <- "loc"
      } else if (method %in% c("join.by.ind", "sidebyside")){
        if(verbose >= 3){
          cat(report(" Joining two genlight datasets with the same individuals but different loci\n"))
        }
        flag <- "ind"
      } else {
        stop(error("Fatal Error: method parameter is deprecated, no longer required. Please remove from function call\n"))
      }
      # Validate the requested join against the data, as the auto-detect
      # path does, so a mismatch fails clearly here rather than
      # cryptically in cbind/rbind
      if (flag == "loc" && !identical(locNames(x1), locNames(x2))) {
        stop(error("Fatal Error: method requests a join by shared loci, but the loci in the two genlight objects do not match\n"))
      }
      if (flag == "ind" && !identical(indNames(x1), indNames(x2))) {
        stop(error("Fatal Error: method requests a join by shared individuals, but the individuals in the two genlight objects do not match\n"))
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
      stop(error("Fatal Error: Individuals or loci in the two files do not match\n"))
    }
}

  # Combine the locus metrics flags: set to 1 only if 1 in both objects;
  # multiply only flags present in both (SNP-style flag data.frames lack
  # OneRatio/PIC and an unconditional product crashes the assignment)
  combine.flags <- function(x, x1, x2, verbose) {
    known.flags <- c("AvgPIC", "OneRatioRef", "OneRatioSnp", "PICRef",
                     "PICSnp", "CallRate", "maf", "FreqHets", "FreqHomRef",
                     "FreqHomSnp", "monomorphs", "OneRatio", "PIC")
    f1 <- x1@other$loc.metrics.flags
    f2 <- x2@other$loc.metrics.flags
    if (!is.null(f1) && !is.null(f2)) {
      for (fl in intersect(known.flags, intersect(names(f1), names(f2)))) {
        x@other$loc.metrics.flags[[fl]] <- f1[[fl]] * f2[[fl]]
      }
    } else {
      if (verbose >= 2) {
        cat(warn(
          "  Warning: Input genlight objects lack locus metrics flags. Flags set to zero\n"
        ))
      }
      for (fl in known.flags) {
        x@other$loc.metrics.flags[[fl]] <- 0
      }
    }
    return(x)
  }
# DO THE JOB --------------
    
if (verbose >= 2) {
  if(flag == "ind"){
    cat(
      report(
        "  Concatenating two genlight objects,",
        x1.name,
        "and",
        x2.name,
        "with the same individuals, different loci\n"
      )
    )
  }
  
  if(flag == "loc"){
    cat(
      report(
        "  Concatenating two genlight objects,",
        x1.name,
        "and",
        x2.name,
        "with the same loci, different individuals\n"
      )
    )
  }
}

if(flag == "ind"){
  if (verbose >= 3) {
    cat("    Number of individuals:", nInd(x1), "\n")
    cat("    First genlight object",
        x1.name,
        "has",
        nLoc(x1),
        "loci\n")
    cat("    Second genlight object",
        x2.name,
        "has",
        nLoc(x2),
        "loci\n")
  }
  
  # Join the two genlight objects (cbind.dartR handles loc.metrics merging)
  x <- cbind(x1, x2)

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
    if (verbose >= 2) {
      cat(
        warn(
          "  Warning: Input genlight objects lack individual metrics\n"
        )
      )
    }
  }

  # Add the loc metrics flags, set to 1 only if 1 in both genlight objects
  if (verbose >= 2) {
    cat(report("  Setting the locus metrics flags\n"))
  }
  x <- combine.flags(x, x1, x2, verbose)
}

if(flag=="loc"){
  if (verbose >= 3) {
    cat("    Number of loci:", nLoc(x1), "\n")
    cat("    First genlight object",
        x1.name,
        "has",
        nInd(x1),
        "individuals\n")
    cat("    Second genlight object",
        x2.name,
        "has",
        nInd(x2),
        "individuals\n")
  }
  
  # Join the two genlight objects
  x <- rbind(x1, x2)

  # Add the individual metrics, combining those of the two objects (rbind
  # does not carry them across)
  if (verbose >= 2) {
    cat(report("  Adding the individual metrics\n"))
  }
  if (!is.null(x1@other$ind.metrics) && !is.null(x2@other$ind.metrics)) {
    x@other$ind.metrics <- plyr::rbind.fill(
      as.data.frame(x1@other$ind.metrics),
      as.data.frame(x2@other$ind.metrics))
  } else {
    if (verbose >= 2) {
      cat(
        warn(
          "  Warning: Input genlight objects lack individual metrics\n"
        )
      )
    }
  }

  # Cater for some individual names being the same in both genlight objects
  if(length(unique(indNames(x))) < length(indNames(x))){
    indNames(x) <- make.unique(indNames(x))
    if (!is.null(x@other$ind.metrics)) {
      x@other$ind.metrics$id <- indNames(x)
    }
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
    if (verbose >= 2) {
      cat(
        warn(
          "  Warning: Input genlight objects lack locus metrics\n"
        )
      )
    }
  }
  
  # Add the loc metrics flags, set to 1 only if 1 in both genlight objects
  if (verbose >= 2) {
    cat(report("  Setting the locus metrics flags\n"))
  }
  x <- combine.flags(x, x1, x2, verbose)
}

    if (verbose >= 3) {
        cat("    Combined genlight object has", nInd(x), "individuals\n")
        cat("    Combined genlight object has", nLoc(x), "loci\n")
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
  
  # FLAG SCRIPT END ---------------
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  # End block -----------------------
  
  return(x)
}
