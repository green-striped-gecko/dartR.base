#' @name gl2snapper
#' @title Converts a genlight object to nexus format suitable for phylogenetic analysis
#'  by Snapper (via BEAUti)
#'  @family linker

#' @description
#' Produces a nexus file contains the SNP calls and relevant PAUP command lines
#' suitable for for the software package BEAUti.

#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile Name of the output file (including extension)
#' [default "snapper.nex"].
#' @param nloc Number of loci to subsample to bring down computational time [default NULL]
#' @param rm.autapomorphies Prune the loci by removing autapomorphies 
#' (not phylogentically informative), that is, SNP polymorphisms limited
#' to only one population [default TRUE].
#' @param outpath Path where to save the output file [default global working 
#' directory or if not specified, tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' 
#' @details {
#'Snapper is a phylogenetic approach implemented in Beauti/BEAST that can
#'handle larger SNP datasets than can SNAPP. This script produces a nexus file 
#'for Beauti which allows options to be set and creates the xml file for BEAST.
#'
#'Although improved over SNAPP in terms of computational efficiency,
#'Snapper remains constrained by computational times, even when implemented on high
#'performance computers with parallel processing. Computation time of Snapper 
#'is not particularly sensitive to the number of individuals in each terminal 
#'taxon (= population), but is impacted by the number of populations and the 
#'numbers of loci. 
#'
#'Particular attention needs to be directed at amalgamating populations
#'where their joint taxonomic identity is without question and reducing the 
#'number of SNP loci by prudent filtering. Removing monomorphic loci, increasing
#'the reliability of loci (read depth, reproducibility), minimizing missing data
#'(call rate), removing multiple SNPs in a single sequence tag (secondaries) 
#'are all good options. You may wish also to remove SNP loci that are polymorphic
#'within only a single population (autapomorphies) are all options for reducing the 
#'number of loci (rm.autapomorphies=TRUE)
#'
#'If computational time is still an issue (say requiring one month for a single run),
#'following strategy is recommended. First, subsample the loci to say 100 and test the process to ensure
#'there are no syntax or other issues (nloc=100). You do not want to wait several days or weeks running
#'the full dataset to discover a simple syntax error or incorrectly specified parameter. 
#'Second consider running BEAST on a platform that allows multi-threading as this will
#'dramatically reduce compute time. Note that adding threads does not always improve
#'computational time. The optimal number of threads depends on the particular analysis.
#'This means you have to experiment to find out how many threads give the best
#'performance foa a computer cycle.
#'
#'Also, there is an overheads cost in using many threads. Instead, you could run independent 
#'snapper analyses and combine resulting log and tree files. For example, if the optimal number of threads is 8 
#' (adding more threads reduces speed), but 8 threads gives marginal improvement over 4 threads, you can run 2 
#'chains with 4 threads each instead and (after getting through burn-in) then combine results and get a better result than running a 
#'single chain at 8 threads. A bonus benefit from running multiple chains is that you can verify that the MCMC ends up with the same 
#'posterior distribution each time.
#'
#'If computational time is still an issue, run the #'analysis on a series of subsamples of loci 
#'(say nloc=300) to see if a consistent topology is obtained,
#'then adopt that topology as the final result.
#'
#'Note that there is a cost to manipulating your data to achieve reasonable
#'computation times. Omission of sequence tags that are invariant during the SNP calling 
#'process, removal of monomorphic loci generated during taxon selection, and removing
#'autapomorphic loci will all affect branch lengths, perhaps differentially, and so 
#'compromise branch lengths and estimates of divergence times. Fortunately, the topology should be little affected.
#'
#'Finally, the analysis relies on certain assumptions. First is that the structure is one of a
#'bifurcating tree and not a network. One needs to assign individuals to populations in advance
#'of the analysis, confident that they are discrete entities and free of horizontal 
#'transfer. A second assumption is that the loci scored for SNPs are assorting 
#'independently. This is probably a reasonable assumption for
#'SNPs derived from sparse representational sampling (e.g. DArT), but if dense SNP arrays are being
#'used, then some form of thinning will be required. Of course, multiple SNPs on a single
#'sequence tag will be linked, so filtering all but one SNP per sequence tag is required
#'(gl.filter.secondaries).
#'
#'The workflow is
#'
#'  \itemize{
#'  \item "1" Execute gl2snapper()
#'  \item "2" Install beast2
#'  \item "3" Run beauti in the beast2 bin
#'  \item "4" Set the template to snapper [File | Template | Snapper]
#'  \item "5" Load the nexus file produced by gl2snapper()
#'  \item "6" Select and set the parameters you consider appropriate
#'  \item "7" Save the xml file [File | Save As]
#'  \item "8" Run beast
#'  \item "9" Load the xml file and execute
#'  \item "10" When beast is finished, examine the diagnostics with Tracer
#'  \item "11" Visualize the resultant trees using DensiTree and FigTree. 
#'  }
#'  
#'  If using the command line to run beast, the command is beast -threads myxmlfile. Progress can be monitored
#'  with awk (awk '\{print $1, $2\}' snapper.log |tail). When beast is finished, transfer the log and tree
#'  files to a windows platform and use Tracer, DensiTree and FigTree as above.
#'
#'gl2snapper does not work with SilicoDArT data.
#'}
#' 
#' @author Custodian: Arthur Georges (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#' x <- gl.filter.monomorphs(testset.gl)
#' gl2snapper(x, outfile="test.nex", outpath=tempdir())
#' 
#' @references Bryant, D., Bouckaert, R., Felsenstein, J., Rosenberg, N.A. and
#' RoyChoudhury, A. (2012). Inferring species trees directly from biallelic
#' genetic markers: bypassing gene trees in a full coalescent analysis.
#'  Molecular Biology and Evolution 29:1917-1932.
#'  
#'  Rambaut A, Drummond AJ, Xie D, Baele G and Suchard MA (2018) Posterior 
#'  summarisation in Bayesian phylogenetics using Tracer 1.7. Systematic Biology. 
#'  syy032. doi:10.1093/sysbio/syy032
#'  
#' @import stats
#' @export
#' @return  returns no value (i.e. NULL)

gl2snapper <- function(x,
                     outfile = "snapper.nex",
                     rm.autapomorphies=FALSE,
                     nloc=NULL,
                     outpath = NULL,
                     verbose = NULL) {

    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # SET WORKING DIRECTORY
    outpath <- gl.check.wd(outpath,verbose=0)
    outfilespec <- file.path(outpath, outfile)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, accept=c("genlight","SNP"),verbose = verbose)
    
    # ERROR CHECKING
    # Spaces in individual names
      if(any(grepl(" ", indNames(x)))){
        cat(warn(" Warning: Names for some entities contain spaces, converted to underscore\n"))
        indNames(x) <- gsub(" ", "_", indNames(x))
      }
    # Duplicate individual names
      if(length(indNames(x))!=length(unique(indNames(x)))){
        cat(warn(" Warning: Names for some entities are not unique, rendering unique\n"))
        dupes <- duplicated(indNames(x)) | duplicated(indNames(x), fromLast = TRUE)
        if(any(dupes)) {
          suffixes <- stats::ave(indNames(x)[dupes], indNames(x)[dupes], FUN = seq_along)
          indNames(x)[dupes] <- paste0(indNames(x)[dupes], "_", suffixes)
        }
      }
    
    # PREPROCESSING
    # Remove autapomorphies
    if(rm.autapomorphies){
      freqs <- gl.allele.freq(x,by = "popxloc")
      subset <- freqs[freqs$frequency<0.99,]
      subset <- subset[subset$frequency>0.0001,]
      tbl <- table(subset$locus)
      single_count_loci <- names(tbl[tbl == 1])
      x <- gl.drop.loc(x,loc.list = single_count_loci)
      cat(report(" Autapomorphic loci removed:",length(single_count_loci),"\n"))
    }
    
    # Subsample loci 
    if(!is.null(nloc)){
      x <- gl.subsample.loc(x,n=nloc,replace=FALSE,verbose=0)
      cat(report(" Loci subsampled at random:",nloc,"retained\n"))
    }

    # DO THE JOB
    
    if (verbose >= 2) {
        cat(paste(
            report(
                "  Extracting SNP data and creating records for each individual\n"
            )
        ))
    }
    
    # Extract the reference base and the alternate base for each locus (excuse the contortion)
    m <- as.matrix(x)
    m[is.na(m)] <- "?"
    colnames(m) <- NULL
    df <- data.frame(m)
    df <- cbind(indNames(x), pop(x), df)
    indlabels <- df[, 1]
    poplabels <- df[, 2]
    
    # Create the snapp file
    if (verbose > 1) {
        cat(report(
            paste("  Writing results to nexus file", outfilespec, "\n")
        ))
    }
    
    sink(outfilespec)
    
    cat("#nexus\n")
    cat("BEGIN DATA;\n")
    cat(paste0("     dimensions ntax = ", nInd(x), " nchar = ", nLoc(x), " ;\n"))
    cat("     format datatype=integerdata missing=? symbols=\"012\";\n")
    cat("matrix\n")
    for (i in 1:nInd(x)) {
        cat(paste0(poplabels[i], "_", indlabels[i]))
        cat("  ")
        cat(m[i, ], sep = "")
        cat("\n")
    }
    cat(";\n")
    cat("end;\n\n")
    
    sink()
    
    if (verbose > 2) {
        cat(paste(
            report("    Records written to", outfilespec, ":", nInd(x), "\n")
        ))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(NULL)
    
}
