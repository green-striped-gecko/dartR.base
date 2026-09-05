#' @name gl2paup.svdquartets
#' @title Converts a genlight object to nexus format PAUP SVDquartets
#' @family linkers

#' @description
#' The output nexus file contains the SNP data in one of two forms, depending
#' upon what you regard as most appropriate. One form, that used by Chifman and
#'  Kubatko, has two lines per individual, one providing the reference SNP the
#'   second providing the alternate SNP (method=1).

#'   A second form, recommended by Dave Swofford, has a single line per
#'   individual, resolving heterozygous SNPs by replacing them with standard
#' ambiguity codes (method=2).

#' If the data are tag presence/absence, then method=2 is assumed and the
#' nexus file is written with datatype = standard (0/1 characters).

#' @param x Name of the genlight object containing the SNP data or tag P/A data
#' [required].
#' @param outfile File name of the output file (including extension)
#' [default 'svd.nex']. The PAUP log and tree file names inside the file's
#' paup block are derived from this name.
#' @param outpath Path where to save the output file [default global working
#' directory or if not specified, tempdir()].
#' @param method Method = 1, nexus file with two lines per individual; method =
#'  2, nexus file with one line per individual, ambiguity codes for SNP
#'  genotypes, 0 or 1 for presence/absence data [default 2].
#' @param nbootstraps Number of bootstrap replicates [default 10000]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#' @return returns no value (i.e. NULL)
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @references Chifman, J. and L. Kubatko. 2014. Quartet inference from SNP data
#' under the coalescent. Bioinformatics 30: 3317-3324
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' gg <- testset.gl[1:20,1:100]
#' gg@other$loc.metrics <- gg@other$loc.metrics[1:100,]
#' gl2paup.svdquartets(gg, outpath=tempdir(),nbootstraps=100)
#'
#' @export

gl2paup.svdquartets <- function(x,
                           outfile = "svd.nex",
                           outpath = NULL,
                           method = 2,
                           nbootstraps=10000,
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
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    if (!is(x, "dartR")) {
      class(x) <- "dartR"  
      if (verbose>2) {
        cat(warn("Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"))
        cat(warn("                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR, 
                 please use function dartR2gl\n"))
      }
    }
    
    # Check for monomorphic loci
    # Tolerate objects not built by dartR, whose loc.metrics.flags are absent (F3)
    monomorphs.flag <- x@other$loc.metrics.flags$monomorphs
    if (is.null(monomorphs.flag) || isFALSE(monomorphs.flag)) {
        if (verbose >= 2) {
            cat(warn(
                "  Warning: genlight object may contain monomorphic loci\n"
            ))
        }
    }

    # FUNCTION SPECIFIC ERROR CHECKING

    # method = 0 was silently accepted by the old range check (F5)
    if (!method %in% c(1, 2)) {
        # Results-affecting coercion: warn at verbose >= 1 (VRB4)
        if (verbose >= 1) {
            cat(
                warn(
                    "  Warning: method must be either 1 or 2, set to 2, one line per individual using ambiguity codes for SNP data, 0 or 1 for tag P/A data \n"
                )
            )
        }
        method <- 2
    }

    # The roxygen promises "If the data are tag presence/absence, then
    # method=2 is assumed"; without this coercion method=1 crashed on the
    # unbuilt refseq/altseq (documented assumption, implemented with F5)
    if (datatype == "SilicoDArT" && method != 2) {
        if (verbose >= 1) {
            cat(warn(
                "  Warning: method=2 is assumed for tag presence/absence data\n"
            ))
        }
        method <- 2
    }

    # Render lables consistent with PAUP
    pop(x) <- gsub(" ", "_", pop(x))
    pop(x) <- gsub("\\(", "_", pop(x))
    pop(x) <- gsub(")", "_", pop(x))
    # Individual names get the same substitution set as populations (F7)
    indNames(x) <- gsub(" ","_",indNames(x))
    indNames(x) <- gsub("\\(", "_", indNames(x))
    indNames(x) <- gsub(")", "_", indNames(x))

    # Sort on population once, ahead of the ploidy branch, so SilicoDArT
    # taxpartition ranges correspond to the matrix rows (F2, F4)
    x <- gl.sort(x, sort.by = "pop", verbose = 0)

    # DO THE JOB
    
    ######## SNP data
    
    if (all(x@ploidy == 2)) {
        if (verbose >= 2) {
            cat(report(
                paste(
                    "    Extacting SNP bases and creating records for each individual\n"
                )
            ))
        }
        
        # Extract the reference base and the alternate base for each locus snp <- as.character(x@other$loc.metrics$SNP) snp <-
        # strsplit(snp,':') snp <- unlist(lapply(snp, `[`, 2)) snp <- strsplit(snp,'>')
        ref <-
            unname(sapply(x@loc.all, function(x)
                strsplit(x, split = "/")[1][[1]][1]))
        # ref <- unlist(lapply(snp, `[`, 1))
        alt <-
            unname(sapply(x@loc.all, function(x)
                strsplit(x, split = "/")[1][[1]][2]))
        # alt <- unlist(lapply(snp, `[`, 2))
        
        df <- data.frame(as.matrix(x))
        df <- cbind(indNames(x), pop(x), df)
        df <- df[order(df$pop), ]
        indlabels <- df[, 1]
        poplabels <- df[, 2]
        m <- df[, 3:(nLoc(x) + 2)]
        
        # Create a lookup table for the ambiguity codes
        if (verbose >= 2) {
            cat(report(
                paste("    Creating lookup table for ambiguity codes\n")
            ))
        }
        # A T G C A A W R M) T W T K Y G R K G S C M Y S C
        
        conversion <-
            matrix(
                c(
                    "A",
                    "W",
                    "R",
                    "M",
                    "W",
                    "T",
                    "K",
                    "Y",
                    "R",
                    "K",
                    "G",
                    "S",
                    "M",
                    "Y",
                    "S",
                    "C"
                ),
                nrow = 4,
                ncol = 4
            )
        colnames(conversion) <- c("A", "T", "G", "C")
        rownames(conversion) <- colnames(conversion)
        
        # Add individual labels to sequences
        refseq <- array(data = NA, dim = nInd(x))
        altseq <- array(data = NA, dim = nInd(x))
        ambseq <- array(data = NA, dim = nInd(x))
        if (method == 1) {
            for (ind in 1:nInd(x)) {
                refseq[ind] <- paste0(indlabels[ind], "_1   ")
                altseq[ind] <- paste0(indlabels[ind], "_2   ")
            }
        } else {
            for (ind in 1:nInd(x)) {
                ambseq[ind] <- paste0(indlabels[ind], "   ")
            }
        }
        
        # progressively add the bases
        if (verbose >= 2) {
            cat(report(
                paste(
                    "    Applying ambiguity codes (if method=2) and reconstructing sequences\n"
                )
            ))
        }
        
        for (ind in 1:nInd(x)) {
            for (loc in 1:nLoc(x)) {
                if (is.na(m[ind, loc])) {
                    refseq[ind] <- paste0(refseq[ind], "?")
                    altseq[ind] <- paste0(altseq[ind], "?")
                    ambseq[ind] <- paste0(ambseq[ind], "?")
                } else if (m[ind, loc] == 0) {
                    refseq[ind] <- paste0(refseq[ind], ref[loc])
                    altseq[ind] <- paste0(altseq[ind], ref[loc])
                    ambseq[ind] <- paste0(ambseq[ind], ref[loc])
                } else if (m[ind, loc] == 1) {
                    refseq[ind] <- paste0(refseq[ind], ref[loc])
                    altseq[ind] <- paste0(altseq[ind], alt[loc])
                    ambseq[ind] <-
                        paste0(ambseq[ind], conversion[ref[loc], alt[loc]])
                } else if (m[ind, loc] == 2) {
                    refseq[ind] <- paste0(refseq[ind], alt[loc])
                    altseq[ind] <- paste0(altseq[ind], alt[loc])
                    ambseq[ind] <- paste0(ambseq[ind], alt[loc])
                } else {
                    cat(error(
                        "Fatal Error: SNP must be scored as 0, 1, 2 or NA\n"
                    ))
                    stop()
                }
            }
        }
    }
    
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
    # Guard for a single population: 2:length(b) would run 2:1 and index b[0] (F4)
    if (length(b) > 1) {
        for (i in 2:length(b)) {
            b[i] <-b[i] + b[i - 1]
            a[i] <-b[i - 1] + 1
        }
    }
    plabels <- unique(poplabels)

    # Derive the PAUP log/tree filenames from outfile (F6)
    fileprefix <- sub("\\.[^.]*$", "", outfile)
    
    # Create the svd file
    if (verbose > 1) {
        cat(report(
            paste("    Writing results to svd nexus file", outfilespec, "\n")
        ))
    }
    
    sink(outfilespec)
    cat("#nexus\n")
    cat("BEGIN DATA;\n")
    if (method == 1) {
        cat(paste0("     dimensions ntax = ", 2 * nInd(x), " nchar = ", nLoc(x), " ;\n"))
    } else {
        cat(paste0("     dimensions ntax = ", nInd(x), " nchar = ", nLoc(x), " ;\n"))
    }
    # SilicoDArT matrices hold 0/1 characters, outside the dna alphabet (F1)
    if (datatype == "SilicoDArT") {
        cat("     format datatype = standard symbols = \"01\" gap = - ;\n\n")
    } else {
        cat("     format datatype = dna gap = - ;\n\n")
    }
    cat("matrix\n")
    if (method == 1) {
        for (i in 1:nInd(x)) {
            cat(paste0(refseq[i], "\n"))
            cat(paste0(altseq[i], "\n"))
        }
    } else {
        for (i in 1:nInd(x)) {
            cat(paste0(ambseq[i], "\n"))
        }
    }
    cat(";\n")
    cat("end;\n\n")
    cat("begin sets;\n")
    cat("    taxpartition pops =\n")
    # With a single population only the terminal line is written (F4)
    if (length(plabels) > 1) {
        for (i in 1:(length(plabels) - 1)) {
            cat(paste0("        ", plabels[i], " : ", a[i], "-", b[i], ",\n"))
        }
    }
    cat("       ", paste0(plabels[length(plabels)], " : ", a[length(plabels)], "-", b[length(plabels)], ";\n"))
    cat("end;\n\n")
    cat("begin paup;\n")
    cat(paste0("log file=",fileprefix,".txt;\n"))
    cat("lset nthreads=3;\n")
    cat("SVDQuartets\n")
    cat("    evalQuartets=random\n")
    cat("    speciesTree=yes\n")
    cat("    partition=pops\n")
    cat("    bootstrap=standard\n")
    cat("    nreps=",nbootstraps,"\n")
    cat("    ambigs=distribute\n")
    cat(paste0("    treeFile=",fileprefix,".tre;\n"))
    cat(paste0("savetrees file=",fileprefix,"_boot.tre from=1 to=1 maxdecimals=2;\n"))
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

    return(invisible(NULL))

}
