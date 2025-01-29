#' @name gl2parsimony
#' @title Converts a genlight object to nexus format for parsimony phylogeny
#' analysis in PAUP.
#' @family linkers

#' @description
#' The output nexus file contains the SilicoDArT data as a single line per
#' individual wrapped in the appropriate nexus commands. Pop Labels are
#' used to define taxon partitions.

#' @param x Name of the genlight object containing the SilicoDArT data
#' [required].
#' @param outfile File name of the output file (including extension)
#' [default 'parsimony.nex'].
#' @param outpath Path where to save the output file [default global working 
#' directory or if not specified, tempdir()].
#' @param nreps Number of bootstrap replicates [default 1000]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity]
#' 
#' @author Custodian: Arthur Georges (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#' gg <- testset.gs[1:20,1:100]
#' gg@other$loc.metrics <- gg@other$loc.metrics[1:100,]
#' gl2parsimony(gg, outpath=tempdir(),nreps=100)
#' 
#' @export
#' @return  returns no value (i.e. NULL)

gl2parsimony <- function(x,
                           outfile = "parsimony.nex",
                           outpath = NULL,
                           nreps=1000,
                           verbose = NULL) {
   
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # SET WORKING DIRECTORY
    outpath <- gl.check.wd(outpath,verbose=0)
    outfilespec <- file.path(outpath, outfile)
     
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, accept="SilicoDArT",verbose = verbose)

    # Check for monomorphic loci
    
    if (!x@other$loc.metrics.flags$monomorphs) {
        cat(warn(
            "  Warning: genlight object may contain monomorphic loci\n"
        ))
    }
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    # Render lables consistent with PAUP
    pop(x) <- gsub(" ", "_", pop(x))
    pop(x) <- gsub("\\(", "_", pop(x))
    pop(x) <- gsub(")", "_", pop(x))
    indNames(x) <- gsub(" ","_",indNames(x))
    
    # DO THE JOB
    
    ######## SilicoDArT data
    
        if (verbose >= 2) {
            cat(report(
                paste(
                    "    Extacting presence-absence data and creating records for each individual\n"
                )
            ))
        }

      # Sort the data on population
        if (verbose >= 2) {
            cat(report(paste("    Sorting ....\n")))
        }
        df <- data.frame(as.matrix(x))
        df <- cbind(indNames(x), pop(x), df)
        df <- df[order(df$pop), ]
        indlabels <- df[, 1]
        poplabels <- df[, 2]
        m <- df[, 3:(nLoc(x) + 2)]
        
    if (all(x@ploidy == 1)) {
        # progressively add the scores (0 0r 1 or NA)
        if (verbose >= 2) {
            cat(report(
                paste("  Constructing genotypes for each individual\n")
            ))
        }
        
        str <- array(NA, nInd(x))
        for (ind in 1:nInd(x)) {
            str[ind] <-
                paste(as.character(as.matrix(x)[ind, ]),
                      collapse = "",
                      sep = "")
            str[ind] <- gsub("NA", "?", str[ind])
            str[ind] <- paste(indNames(x)[ind], "   ", str[ind])
        }
        ambseq <- str
        poplabels <- pop(x)
    }
    
    # Create the taxpartition (popname : 25-60)
    if (verbose >= 2) {
        cat(report(paste("    Creating partition table\n")))
    }
    a <- array(data = NA, dim = length(poplabels))
    b <- array(data = NA, dim = length(poplabels))
    a[1] <- 1
    b <- table(poplabels)
    for (i in 2:length(b)) {
        b[i] <-b[i] + b[i - 1]
        a[i] <-b[i - 1] + 1
    }
    plabels <- unique(poplabels)
    
    # Create the parsimony file
    if (verbose > 1) {
        cat(report(
            paste("    Writing results to parsimony nexus file", outfilespec, "\n")
        ))
    }
    
    sink(outfilespec)
    cat("#nexus\n")
    cat("BEGIN DATA;\n")
    cat(paste0("     dimensions ntax = ", nInd(x), " nchar = ", nLoc(x), " ;\n"))
    cat("     format datatype = dna gap = - ;\n\n")
    cat("matrix\n")
    for (i in 1:nInd(x)) {
      cat(paste0(ambseq[i], "\n"))
    }
    cat(";\n")
    cat("end;\n\n")
    cat("begin sets;\n")
    cat("    taxpartition pops =\n")
    for (i in 1:(length(plabels) - 1)) {
        cat(paste0("        ", plabels[i], " : ", a[i], "-", b[i], ",\n"))
    }
    cat("       ", paste0(plabels[length(plabels)], " : ", a[length(plabels)], "-", b[length(plabels)], ";\n"))
    cat("end;\n\n")
    cat("begin paup;\n")
    cat("log file=parsimony.txt;\n")
    cat("lset nthreads=3;\n")
    cat("set criterion=parsimony;\n")
    cat("hsearch addseq=random nreps=100 swap=tbr multrees=yes;\n")
    cat("bootstrap nreps=",nreps," search=heuristic addseq=random;\n")
    cat("savetrees file=parsimony_boot.tre from=1 to=1 maxdecimals=2;\n")
    cat("log stop;\n")
    cat("quit;\n")
    cat("end;\n")

    sink()
    
    if (verbose > 2) {
        cat(report(
            paste("    Records written to", outfilespec, ":", nInd(x), "\n")
        ))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(NULL)
    
}
