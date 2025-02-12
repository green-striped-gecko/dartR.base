#' @name gl2paup.parsimony
#' @title Converts a genlight object to nexus format for parsimony phylogeny
#' analysis in PAUP and, optionally produces accompanying files for parallel processing.
#' @family linkers

#' @description
#' The output nexus file contains the SilicoDArT data as a single line per
#' individual wrapped in the appropriate nexus commands. Pop Labels are
#' used to define taxon partitions.
#' 
#' If out.type="bash", the function produces a series of files in support of an
#' analysis taking advantage of multi-threading and parallel processing.

#' @param x Name of the genlight object containing the SilicoDArT data
#' [required].
#' @param outfileprefix A prefix to use for file names of the output files
#' [default 'parsimony'].
#' @param outpath Path where to save the output file [default global working 
#' directory or if not specified, tempdir()].
#' @param out.type Specify the type of output file. Can be 'standard' (consensus tree)
#' or 'newick' (newick) or 'bash' [default 'standard']
#' @param tip.labels Specify whether the terminals should be labelled with the
#' individual labels ('ind'), the population labels ('pop') or both ('indpop') 
#' [default 'ind']
#' @param nreps Specify the number of replicate analyses to run in search of
#' the shortest tree [default 100]
#' @param nbootstraps Number of bootstrap replicates [default 1000]
#' @param ncpus Number of cores to use for parallel processing [default 1]
#' @param mem Memory to use for each process [default 4Gb per core]
#' @param server If out.type='bash', provide the name of the linux server [default 'gadi']
#' @param base.dir.name Name of the base directory on the server to act as the workspace [default NULL]
#' @param test If TRUE, the analysis will run with a small subset of the data [default FALSE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity]
#' 
#' @details
#' Additional details: This script only applies to SilicoDArT data. The output file
#' is the name of the file PAUP will use to deliver the results of the analysis, in
#' the directory specified by outpath. 
#' 
#' The output type (out.type) can be 'standard' which uses
#' default PAUP parameters to construct the boot.tre file. Or it can be 'newick' to add the 
#' parameter format=newick whereby the boot.tre file contains the final tree in newick
#' format. This is useful for passing the results to a tree graphics program such as
#' Mega 11 to format the tree for publication. Or it can be 'bash' which creates
#' a number of files to facilitate parallel processing on a supercomputer.
#' 
#' The parameter nreps specifies the number of replicates to run in search of the
#' shortest tree in each bootstrap iteration. The default is 100.
#' 
#' The parameter nbootstraps specifies tne number of bootstrap replicates to run to generate 
#' a measure of node support. The default is 1000. The companion parameter ncpus specifies how
#' many cpus to use for parallel processing when out.type='bash'. The default is 1. 
#' Note that the number of cpus must divide evenly into the number of bootstrap replicates.
#' 
#' The parameter tip.labels specifies whether the terminals in the tree should be
#' labelled with the individual names, or the population names (multiple tips will
#' have the same label -- which can cause problems at the point of generating a 
#' consensus tree), or a combination of the two. Including the population name
#' in the terminal tip labels will assist in collapsing the tree to have populations
#' as the terminals after checking fidelity of populations to supported clades. This
#' can be done in Mega 11.
#' 
#' The parameter 'server' is to allow for future development as users modify the bash
#' scripts to suit other multitasking environments. This script works only for the
#' Gadi server on the Australian National Computing Infrastructure (NCI).
#' 
#' If test=TRUE, the data will be subsetted heavily on numbers of loci, numbers individuals,
#' bootstrap replicates and number of replicates for branch swapping. This is used to test
#' the job run without expenditure of the resources required for the full job.
#' 
#' 
#' @author Custodian: Arthur Georges (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#' gg <- testset.gs[1:20,1:100]
#' gg@other$loc.metrics <- gg@other$loc.metrics[1:100,]
#' gl2paup.parsimony(gg,outfile="test.nex",outpath=tempdir(),nreps=1,nbootstraps=10)
#' gl2paup.parsimony(gg,outfile="test.nex",out.type="newick",outpath=tempdir(),nreps=1,nbootstraps=10)
#' 
#' @export
#' @return  returns no value (i.e. NULL)

gl2paup.parsimony <- function(x,
                         outfileprefix = "parsimony",
                         outpath = NULL,
                         out.type="standard",
                         tip.labels="ind",
                         nreps=100,
                         nbootstraps=1000,
                         ncpus=1,
                         mem=4,
                         server="gadi",
                         base.dir.name=NULL,
                         test=FALSE,
                         verbose = NULL) {
   
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # SET WORKING DIRECTORY
    outpath <- gl.check.wd(outpath,verbose=0)
    outfilespec <- file.path(outpath, outfileprefix)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, accept="SilicoDArT",verbose = verbose)
    
    if (!is(x, "dartR")) {
      class(x) <- "dartR"  
      if (verbose>2) {
        cat(warn("Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"))
        cat(warn("                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR, 
                 please use function dartR2gl\n"))
      }
    }
    
    # IF TEST=TRUE, SUBSET THE ANALYSIS
      if(test==TRUE){
        if(verbose>=2){cat(report("  Test run only\n"))}
        nreps <- 1
        nbootstraps <- 10
        ncpus <- 10
        mem <- ncpus*4
        qtl <- quantile((table(pop(x))),0.75)
        qtl.names <- names(table(pop(x))[table(pop(x))>=qtl])
        x <- gl.keep.pop(x, pop.list=qtl.names,mono.rm=TRUE,verbose=0)
        if(nLoc(x)>1000){x <- gl.subsample.loc(x,n=1000,replace=FALSE,verbose=0)}
        x <- gl.subsample.ind(x,n=qtl,replace=FALSE,verbose=0)
        if(verbose>=2){cat(report("  Test run on",nInd(x),"individuals in",nPop(x),"populations scored for",nLoc(x),"SNP loci\n"))}
        if(verbose>=2){cat(report("  PAUP parameters:",nbootstraps,"bootstrap replicates split over",ncpus,"CPUs","with",nreps,"heuristic tree searches in each iteration\n"))}
      }


    # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x,verbose=0)
    if(nLoc(x) != nLoc(tmp)){
      cat(warn(
        "  Warning: genlight object may contain monomorphic loci\n"
      ))
    }
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if(!out.type %in% c("standard","newick","bash")){
      cat(warn("Warning: Output format (out.type) must be one of 'standard' or 'newick' or 'bash', set to 'standard'\n"))
      out.type <- "standard"
    }
    
    if(out.type=='bash'){
      if((nbootstraps %% ncpus) != 0){
        cat(error("Fatal Error: The number of bootstraps",nbootstraps,"must be a multiple of the number of cpus",ncpus,"\n"))
      }
    }
    
    if(out.type=="bash"){
      if(is.null(base.dir.name)){
        cat(error("Fatal Error: a base name for the working directory on",server,"is required\n"))
        stop()
      }
    }
    
    if(tip.labels=="popind"){tip.labels <- "indpop"}
    if(!tip.labels %in% c("ind","pop","indpop")){
      cat(warn("Warning: Tip labeles (tip.labels) must be one of 'ind','pop' or 'indpop', set to 'ind'\n"))
      tip.labels <- "ind"
    }
    
    if(test==FALSE){
      if(nreps < 100){
        cat(warn("Warning: Parameter nreps recommended to be 100 or greater\n"))
      }
      if(nbootstraps < 1000){
        cat(warn("Warning: Parameter nbootstraps recommended to be 1000 or greater\n"))
      }
    }
    
    # Render lables consistent with PAUP
    pop(x) <- as.character(pop(x))
    pop(x) <- gsub(" ", "_", pop(x))
    pop(x) <- gsub("\\(", "_", pop(x))
    pop(x) <- gsub(")", "_", pop(x))
    pop(x) <- sub("_$", "", pop(x))
    indNames(x) <- gsub(" ","_",indNames(x))
    indNames(x) <- gsub("\\(", "_", indNames(x))
    indNames(x) <- gsub(")", "_", indNames(x))
    indNames(x) <- sub("_$", "", indNames(x))
    
    # DO THE JOB
    
    ###############################################################################################
    ##### out.type='standard' or out.type='newick' or out.type='bash'
    ###############################################################################################
    
        if (verbose >= 2) {
            cat(report(
                paste(
                    "    Extacting presence-absence data and creating records for each individual\n"
                )
            ))
        }

    # Sort the data on population
        x <- gl.sort(x,sort.by="pop",verbose=0)

    if (all(x@ploidy == 1)) {
        # progressively add the scores (0 0r 1 or NA)
        if (verbose >= 2) {
            cat(report(
                paste("  Constructing genotypes for each individual ... coffee time\n")
            ))
        }
      
        # First convert indNames
        tmp <- x
        if(tip.labels=="pop"){
          indNames(tmp) <- pop(x)
        } else if(tip.labels=="indpop" | tip.labels=="popind"){
          indNames(tmp) <- paste0(pop(x),"_",indNames(x))
        }
        
        # Add to sequence strings
        str <- array(NA, nInd(tmp))
        for (i in 1:nInd(tmp)) {
            str[i] <-
                paste(as.character(as.matrix(tmp)[i, ]),
                      collapse = "",
                      sep = "")
            str[i] <- gsub("NA", "?", str[i])
            str[i] <- paste(indNames(tmp)[i], "   ", str[i])
        }
        ambseq <- str
        poplabels <- pop(x)
    }
    
    # Create the taxpartition (popname : 25-60)
    if (verbose >= 2) {
        cat(report(paste("    Creating taxpartition table\n")))
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
    outfilespec1 <- paste0(outfilespec,"_bootstrap.nex")

    sink(outfilespec1)
    cat("#NEXUS\n")
    cat("BEGIN DATA;\n")
    cat(paste0("     dimensions ntax = ", nInd(x), " nchar = ", nLoc(x), " ;\n"))
    cat("     format datatype = standard gap = - ;\n\n")
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
    cat(paste0("log file=",outfileprefix,".log;\n"))
    cat("set criterion=parsimony storebrlens=yes increase=auto storetreewts=yes;\n")

    cat("set autoclose=yes;\n")
    #cat("hsearch addseq=random nreps=",nreps," swap=tbr multrees=yes;\n")
    if(out.type=="bash"){
      cat(paste0("bootstrap nreps=",nbootstraps/ncpus," search=heuristic / start=stepwise addseq=random nreps=",nreps," swap=TBR;\n"))
    } else {
      cat(paste0("bootstrap nreps=",nbootstraps," search=heuristic brlens=yes / start=stepwise addseq=random nreps=",nreps," swap=TBR;\n"))
    }
    cat(paste0("savetrees file=",outfileprefix,"_bootstrap.tre;\n"))
    cat("log stop;\n")
    cat("quit;\n")
    cat("end;\n")

    sink(NULL)
    
    if(verbose > 2){
      cat(report("  Base nexus file written drawing upon",nInd(x),"genotypes from",nPop(x),"populations and",nLoc(x),"loci\n"))
    }
    
    ###############################################################################################
    ##### out.type='bash'
    ###############################################################################################
    if(out.type=="bash"){
    # Create the driver for creating bootstrap nexus files bootstrap1.nex to bootstrapN.nex

    # Define the output file name
    outfilespec2 <- paste0("generator_", outfileprefix, "_bootstraps.sh")

    # Ensure `ncpus` is defined
    num_iterations <- as.integer(ncpus)  # Convert to integer to prevent NA issues

    # Open a connection to write in binary mode (ensures Unix `LF` line endings)
    con <- file(outfilespec2, open = "wb")

    # Write the Bash script using `writeLines()`, ensuring one statement per line
    writeLines("#!/bin/bash", con)
    writeLines("", con)

    writeLines("# Locate the _bootstrap.nex file", con)
    writeLines("NEX_FILE=$(ls *_bootstrap.nex 2>/dev/null)", con)

    writeLines("# Ensure exactly one _bootstrap.nex file exists", con)
    writeLines("if [[ $(echo \"$NEX_FILE\" | wc -l) -ne 1 ]]; then", con)
    writeLines("  echo 'Error: There should be exactly one _bootstrap.nex file.'", con)
    writeLines("  exit 1", con)
    writeLines("fi", con)

    writeLines("# Read the template file into a variable", con)
    writeLines("NEX_CONTENT=$(cat \"$NEX_FILE\")", con)

    writeLines("# Calculate the number of bootstrap iterations dynamically", con)
    writeLines(paste0("NUM_ITERATIONS=", num_iterations), con)

    writeLines("# Loop to generate bootstrap files", con)
    writeLines("for i in $(seq 1 $NUM_ITERATIONS); do", con)
    writeLines("  BOOTSTRAP_FILE=\"bootstrap${i}.nex\"", con)

    writeLines("  # Replace log file and tree file names dynamically", con)
    # writeLines("  MODIFIED_NEX=$(echo \"$NEX_CONTENT\" | \
    # sed 's/log file=[^ ;]*/log file=bootstrap${i}.log;/g' | \
    # sed 's/savetrees file=[^ ;]*/savetrees file=bootstrap'\"${i}\"'.tre;/g')", con)
    writeLines("  MODIFIED_NEX=$(echo \"$NEX_CONTENT\" | \
    sed 's/log file=[^ ;]*/log file=bootstrap'\"$i\"'.log replace;/g' | \
    sed 's/savetrees file=[^ ;]*/savetrees file=bootstrap'\"$i\"'.tre;/g')", con)
    writeLines("", con)

    writeLines("  # Write modified Nexus file", con)
    writeLines("  echo \"$MODIFIED_NEX\" > \"$BOOTSTRAP_FILE\"", con)

    writeLines("  echo \"Generated: $BOOTSTRAP_FILE\"", con)
    writeLines("done", con)

    writeLines("echo \"All bootstrap Nexus files have been created successfully!\"", con)

    # Close the file connection
    close(con)

    # Optional verbosity message
    if (verbose > 2) {
      cat(report("  Generator bash script", outfilespec2, "written to generate", ncpus, "nexus files\n"))
    }

    ############################################################################

    # Create the driver for running the ncpus nexus files in parallel to produce bootstrap1.tre
    # to bootstrap20.tre)

    outfilespec3 <- paste0("generator_", outfileprefix, "_maketrees.sh")

    # Open a connection to write in binary mode (forces Unix line endings)
    con <- file(outfilespec3, open = "wb")

    writeLines("#!/bin/bash\n", con)
    writeLines(paste0("NUM_BOOTSTRAPS=",nbootstraps/ncpus,"\n"), con)

    writeLines("for i in $(seq 1 $NUM_BOOTSTRAPS); do", con)
    writeLines("JOB_FILE=\"bootstrap_job${i}.pbs\"\n", con)  # No escaping needed

    # Begin `cat <<EOT` block correctly (no escaping inside EOT)
    writeLines("cat <<EOT > \"$JOB_FILE\"\n", con)
    writeLines("#!/bin/bash\n", con)

    # PBS directives
    writeLines("#PBS -P xl04", con)
    writeLines("#PBS -q normal", con)
    writeLines(paste0("#PBS -l ncpus=",ncpus), con)
    writeLines(paste0("#PBS -l mem=",mem,"GB"), con)
    writeLines("#PBS -l walltime=48:00:00", con)
    writeLines("#PBS -j oe", con)
    writeLines(paste0("#PBS -o ",base.dir.name,"/pbslogs"), con)
    writeLines(paste0("#PBS -N ", outfileprefix, "_bootstrap_job${i}"), con)  # No escaping needed
    writeLines("#PBS -l storage=gdata/xl04+gdata/if89\n", con)

    # Load PAUP module
    writeLines("module load paup\n", con)

    # Navigate to the working directory
    writeLines("cd /g/data/xl04/ag3760/parsimony/\n", con)

    # Run PAUP with correct variable expansion (NO escaping needed inside EOT)
    writeLines("/g/data/if89/apps/paup/4a168/paup -n bootstrap${i}.nex\n", con)

    # Print job completion message
    writeLines("echo \"Bootstrap job ${i} completed successfully\"\n", con)

    writeLines("EOT\n", con)  # Correctly close EOT (no escaping)

    # Submit the job (escaping needed since it's outside EOT)
    writeLines("qsub $JOB_FILE", con)
    writeLines("echo \"Submitted bootstrap job ${i}\"\n", con)
    writeLines("done", con)
    writeLines("echo \"All jobs have been submitted to the queue for execution!\"\n", con)

    # Close the file connection
    close(con)

    if (verbose > 2) {
      cat(report("  Generator bash script", outfilespec3, "written to generate", ncpus, "tree files\n"))
    }

    ###### WRITE THE SCRIPT TO CALCULATE THE CONSENSUS TREE

    outfilespec4 <- paste0("generator_", outfileprefix, "_consensus.sh")

    # Open a connection to write in binary mode (forces Unix line endings)
    con <- file(outfilespec4, open = "wb")

    # PBS job settings
    writeLines("#PBS -P xl04", con)
    writeLines("#PBS -q normal", con)
    writeLines(paste0("#PBS -l ncpus=",ncpus), con)
    writeLines(paste0("#PBS -l mem=",mem,"GB"), con)
    writeLines("#PBS -l walltime=48:00:00", con)
    writeLines("#PBS -j oe", con)
    writeLines(paste0("#PBS -o ",base.dir.name,"/pbslogs"), con)
    writeLines(paste0("#PBS -N ", outfileprefix, "_bootstrap_job${i}"), con)  # No escaping needed
    writeLines("#PBS -l storage=gdata/xl04+gdata/if89\n", con)

    writeLines("module load paup\n", con)

    # Navigate to working directory
    writeLines(paste0("cd ",base.dir.name,"\n"), con)

    # Create PAUP Nexus file dynamically using a Bash loop
    writeLines("echo 'begin paup;' > final_consensus.nex", con)
    writeLines("echo 'log file=consensus.log;' >> final_consensus.nex", con)
    writeLines("echo 'set criterion=parsimony;' >> final_consensus.nex", con)
    writeLines("echo '' >> final_consensus.nex", con)
    writeLines("echo 'set storebrlens=yes;' >> final_consensus.nex", con)
    writeLines("echo 'set increase=auto;' >> final_consensus.nex", con)
    writeLines("echo '[Load all bootstrap trees dynamically]' >> final_consensus.nex\n", con)

    # Loop over existing bootstrap trees and append to PAUP file
    writeLines("for file in bootstrap*.tre; do", con)
    writeLines("  if [[ -f \"$file\" ]]; then", con)
    writeLines("    echo \"gettrees file=$file mode=7;\" >> final_consensus.nex", con)
    writeLines("  fi", con)
    writeLines("done\n", con)

    # Add consensus command to PAUP script
    writeLines("echo '' >> final_consensus.nex", con)
    writeLines("echo '[Calculate majority-rule consensus tree with bootstrap values]' >> final_consensus.nex", con)
    writeLines("echo 'contree all / strict=no majrule=yes percent=50 treefile=final_consensus.nwk format=newick replace;' >> final_consensus.nex", con)
    writeLines("echo '' >> final_consensus.nex", con)
    writeLines("echo 'log stop;' >> final_consensus.nex", con)
    writeLines("echo 'quit;' >> final_consensus.nex", con)
    writeLines("echo 'end;' >> final_consensus.nex\n", con)

    # Run PAUP with the generated consensus Nexus file
    writeLines("/g/data/if89/apps/paup/4a168/paup -n final_consensus.nex\n", con)

    writeLines("echo 'Consensus tree computation completed successfully'\n", con)

    # Close the file connection
    close(con)

    if (verbose > 2) {
      cat(report("  Generator bash script", outfilespec4, "written to generate and submit the PAUP consensus tree job on",server,"\n"))
    }
    }

    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(NULL)
    
}
