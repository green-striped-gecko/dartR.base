#' Concatenates DArT trimmed sequences and outputs a FASTA file
#'
#' Concatenated sequence tags are useful for phylogenetic methods where
#' information on base frequencies and transition and transversion ratios are
#' required (for example, Maximum Likelihood methods). Where relevant,
#' heterozygous loci are resolved before concatenation by either assigning
#'  ambiguity codes or by random allele assignment.
#'
#' Four methods are employed:
#'
#' Method 1 -- heterozygous positions are replaced by the standard ambiguity 
#' codes. The resultant sequence fragments are concatenated across loci to
#' generate a single combined sequence to be used in subsequent ML phylogenetic
#' analyses.
#'
#' Method 2 -- the heterozygous state is resolved by randomly assigning one or
#' the other SNP variant to the individual. The resultant sequence fragments are
#' concatenated across loci to generate a single composite haplotype to be used
#' in subsequent ML phylogenetic analyses.
#'
#' Method 3 -- heterozygous positions are replaced by the standard ambiguity
#' codes. The resultant SNP bases are concatenated across loci to generate a
#' single combined sequence to be used in subsequent MP phylogenetic analyses.
#'
#' Method 4 -- the heterozygous state is resolved by randomly assigning one or
#' the other SNP variant to the individual. The resultant SNP bases are
#' concatenated across loci to generate a single composite haplotype to be used
#' in subsequent MP phylogenetic analyses.
#'
#' Trimmed sequences for which the SNP has been trimmed out, rarely, by adapter
#'  mis-identity are deleted.
#'
#' The script writes out the composite haplotypes for each individual as a
#' fastA file. Requires 'TrimmedSequence' to be among the locus metrics
#' (\code{@other$loc.metrics}) and information of the type of alleles (slot
#' loc.all e.g. 'G/A') and the position of the SNP in slot position of the
#'  ```genlight``` object (see testset.gl@position and testset.gl@loc.all for
#'  how to format these slots.)
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param method One of 1 | 2 | 3 | 4. Type method=0 for a list of options 
#' [method=1].
#' @param outfile Name of the output file (fasta format) ["output.fasta"].
#' @param outpath Path where to save the output file (set to tempdir by default)
#' @param probar If TRUE, a progress bar will be displayed for long loops 
#' [default = TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log; 3, progress and results summary; 5, full report 
#'  [default 2 or as specified using gl.set.verbosity].
#' @return NULL, invisibly. The FASTA file is written to outpath.
#' @export
#' @importFrom utils combn edit flush.console getTxtProgressBar read.csv setTxtProgressBar txtProgressBar write.csv write.table
#' @importFrom graphics axis barplot box image lines text
#' @importFrom methods new
#' @importFrom stats dist nobs optimize pchisq variable.names
#' @import stringr
#' @author Custodian: Luis Mijangos (Post to
#'  \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'  \donttest{
#'  if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' gl <- gl.filter.reproducibility(testset.gl,t=1)
#' gl <- gl.filter.overshoot(gl,verbose=3)
#' gl <- gl.filter.callrate(testset.gl,t=.98)
#' gl <- gl.filter.monomorphs(gl)
#' gl2fasta(gl, method=1, outfile='test.fasta',verbose=3)
#' }

gl2fasta <- function(x,
                     method = 1,
                     outfile = "output.fasta",
                     outpath = tempdir(),
                     probar = FALSE,
                     verbose = NULL) {
    outfilespec <- file.path(outpath, outfile)

    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    # CHECK IF PACKAGES ARE INSTALLED
    pkg <- "seqinr"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      stop(error("Package", pkg, "needed for this function to work.
                    Please install it.\n"))
    }
    
    # Check monomorphs have been removed up to date
    if (!isTRUE(x@other$loc.metrics.flags$monomorphs)) {
        if (verbose >= 2) {
            cat(
                warn(
                    "  Warning: Dataset contains monomorphic loci which will be included in the output fasta file\n"
                )
            )
        }
    }
    
    if (length(x@other$loc.metrics$TrimmedSequence) != nLoc(x)) {
        stop(
            error(
                "Fatal Error: Data must include Trimmed Sequences for each loci in a column called 'TrimmedSequence' in the @other$loc.metrics slot.\n"
            )
        )
    }
    #if (length(x@position) != nLoc(x)) {
    if (length(x@other$loc.metrics$SnpPosition) != nLoc(x)) {
      stop(
        error(
          "Fatal Error: Data must include position information for each loci in x@other$loc.metrics$SnpPosition.\n"
        )
      )
    }
    if (length(x@loc.all) != nLoc(x)) {
      stop(error(
        "Fatal Error: Data must include type of alleles in the @loc.all slot.\n"))
    }

    # Coerce SnpPosition to numeric up front; a factor or character column
    # would otherwise feed substr with factor level codes and silently
    # corrupt the output sequences
    raw_pos <- x@other$loc.metrics$SnpPosition
    num_pos <- suppressWarnings(as.numeric(as.character(raw_pos)))
    if (sum(is.na(num_pos)) > sum(is.na(raw_pos))) {
        stop(error("Fatal Error: @other$loc.metrics$SnpPosition contains
                    non-numeric values.\n"))
    }
    x@other$loc.metrics$SnpPosition <- num_pos
    
    if (method == 1) {
        if(verbose >=2){cat(
            report(
                "  Assigning ambiguity codes to heterozygote SNPs, concatenating trimmed sequence\n"
            )
        )}
    } else if (method == 2) {
      if(verbose >=2){cat(
            report(
                "  Randomly allocating heterozygotes (1) to homozygotic state (0 or 2), concatenating trimmed sequence\n"
            )
        )}
    } else if (method == 3) {
      if(verbose >=2){cat(
            report(
                "  Assigning ambiguity codes to heterozygote SNPs, concatenating SNPs\n"
            )
        )}
    } else if (method == 4) {
      if(verbose >=2){cat(
            report(
                "  Randomly allocating heterozygotes (1) to homozygotic state (0 or 2), concatenating SNPs\n"
            )
        )
      }
    } else {
        if (verbose >= 2) {
            cat(warn("Method not properly specified\n"))
            cat(warn("  Replace score for heterozygotic loci with:\n"))
            cat(warn(
                "  method=1 -- ambiguity codes, concatenate fragments [default]\n"
            ))
            cat(
                warn(
                    "  method=2 -- random assignment to one of the two homogeneous states, concatenate fragments\n"
                )
            )
            cat(warn(
                "  method=3 -- ambiguity codes, concatenate SNPs only\n"
            ))
            cat(
                warn(
                    "  method=4 -- random assignment to one of the two homogeneous states, concatenate SNPs only\n"
                )
            )
        }
        stop(error("Fatal Error: Parameter method out of range.\n"))
    }
    
    # DO THE JOB
    
    if (verbose >= 2) {
      cat(report(
        paste(
          "  Removing loci for which SNP position is outside the length of the trimmed sequences\n"
        )
      ))
    }
    nloc_start <- nLoc(x)

    # Loci that cannot be rendered at all (missing position or sequence)
    trimmed_all <- as.character(x@other$loc.metrics$TrimmedSequence)
    not_renderable <- which(is.na(x@other$loc.metrics$SnpPosition) |
                                is.na(trimmed_all) | nchar(trimmed_all) == 0)
    if (length(not_renderable) > 0) {
        x <- gl.drop.loc(x, loc.list = locNames(x)[not_renderable],
                         verbose = 0)
    }

    # The overshoot pre-filter (output captured: its own reporting is
    # driven by this function's accounting below)
    invisible(capture.output(x <- gl.filter.overshoot(x, verbose = 0)))

    nloc_removed <- nloc_start - nLoc(x)
    if (nloc_removed > 0 && verbose >= 1) {
        cat(warn(
            "  Warning:", nloc_removed, "of", nloc_start,
            "loci removed by the pre-filter (SNP position missing or
                    outside the trimmed sequence)\n"
        ))
    }
    if (nLoc(x) < 2) {
        stop(error(
            "Fatal Error: only", nLoc(x),
            "locus remains after the overshoot pre-filter. The SnpPosition
                    values are inconsistent with the TrimmedSequence lengths
                    - a known hazard in merged or co-analysed datasets.
                    Check @other$loc.metrics$SnpPosition against
                    nchar(@other$loc.metrics$TrimmedSequence).\n"
        ))
    }
    
    # METHOD = AMBIGUITY CODES
    
    if (method == 1 || method == 3) {
        allnames <- locNames(x)
        snp <- as.character(x@loc.all)
        trimmed <- as.character(x@other$loc.metrics$TrimmedSequence)
        snpmatrix <- as.matrix(x)
        
        # Create a lookup table for the ambiguity codes A T G C A A W R M) T W T K Y G R K G S C M Y S C
        
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
        
        # Extract alleles 1 and 2
        #allelepos <- x@position
        allelepos <- x@other$loc.metrics$SnpPosition
        allele1 <- gsub("(.)/(.)", "\\1", snp, perl = T)
        allele2 <- gsub("(.)/(.)", "\\2", snp, perl = T)
        
        # Prepare the output fastA file
        if (verbose >= 2) {
            cat(report("Generating haplotypes ... This may take some time\n"))
        }
        
        sink(outfilespec)
        on.exit(if (sink.number() > 0) sink(), add = TRUE)

        for (i in 1:nInd(x)) {
            seq <- NA
            for (j in 1:nLoc(x)) {
                if (is.na(snpmatrix[i, j])) {
                    code <- "N"
                } else {
                    if (snpmatrix[i, j] == 0) {
                        a1 <-allele1[j]
                        a2 <-allele1[j]
                    }
                    if (snpmatrix[i, j] == 1) {
                        a1 <-allele1[j]
                        a2 <-allele2[j]
                    }
                    if (snpmatrix[i, j] == 2) {
                        a1 <-allele2[j]
                        a2 <-allele2[j]
                    }
                    code <- conversion[a1, a2]
                }
                snppos <- allelepos[j]
                if (method == 1) {
                    if (code != "N") {
                        seq[j] <-
                            paste0(
                                substr(trimmed[j], 1, snppos),
                                code,
                                substr(trimmed[j], snppos + 2,
                                       nchar(trimmed[j]))
                            )
                    } else {
                        seq[j] <-
                            paste(rep("N", nchar(trimmed[j])), collapse = "")
                    }
                } else if (method == 3) {
                    seq[j] <- code
                }
            }
            # Join all the trimmed sequence together into one long 'composite' haplotype
            result <- paste(seq, sep = "", collapse = "")
            # Write the results to file in fastA format
            cat(paste0(">", indNames(x)[i], "_", pop(x)[i], "\n"))
            cat(result, "\n", sep = "")
        }
        
        # Close the output fastA file
        sink()
    }
    
    # METHOD = RANDOM ASSIGNMENT
    
    if (method == 2 || method == 4) {
        # Randomly allocate heterozygotes (1) to homozygote state (0 or 2)
        genmat <- as.matrix(x)
        r <- nrow(genmat)
        c <- ncol(genmat)
        for (i in 1:r) {
            for (j in 1:c) {
                if (!is.na(genmat[i, j]) && genmat[i, j] == 1) {
                    # Score it 0 or 2
                    genmat[i, j] <- (sample(1:2, 1) - 1) * 2
                }
            }
        }
        
        # Prepare the output fastA file
        if (verbose >= 2) {
            cat(report("Generating haplotypes ... This may take some time\n"))
        }
        
        sink(outfilespec)
        on.exit(if (sink.number() > 0) sink(), add = TRUE)

        # For each individual, and for each locus, generate the relevant haplotype
        seq <- rep(" ", c)
        # pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label='Working ....') getTxtProgressBar(pb)
        for (i in 1:r) {
            for (j in 1:c) {
                # Reassign some variables
                trimmed <-
                    as.character(x@other$loc.metrics$TrimmedSequence[j])
                snp <- x@loc.all[j]
                #snpos <- x@position[j]
                snpos <- x@other$loc.metrics$SnpPosition[j]
                # Shift the index for snppos to start from 1 not zero
                snpos <- snpos + 1
                
                # If the score is homozygous for the reference allele
                if (method == 2) {
                    if (!is.na(genmat[i, j]) && genmat[i, j] == 0) {
                        seq[j] <- trimmed
                    }
                } else if (method == 4) {
                    seq[j] <- stringr::str_sub(trimmed,
                                               start = (snpos),
                                               end = (snpos))
                }
                # If the score is homozygous for the alternate allele
                if (!is.na(genmat[i, j]) && genmat[i, j] == 2) {
                    # Split the trimmed into a beginning sequence, the SNP and an end sequences
                    start <- stringr::str_sub(trimmed, end = snpos - 1)
                    snpbase <- stringr::str_sub(trimmed,
                                         start = (snpos),
                                         end = (snpos))
                    end <- stringr::str_sub(trimmed, start = snpos + 1)
                    # Extract the SNP transition bases (e.g. A and T)
                    state1 <-gsub("(.)/(.)", "\\1", snp, perl = T)
                    state2 <-gsub("(.)/(.)", "\\2", snp, perl = T)
                    # Change the SNP state to the alternate
                    if (snpbase == state1) {
                        snpbase <- state2
                    } else {
                        snpbase <- state1
                    }
                    
                    if (method == 2) {
                        # Paste back to form the alternate fragment
                        target <- paste0(start, snpbase, end)
                        # Remove adaptors and save the trimmed alternate sequence
                        seq[j] <-
                            stringr::str_sub(target,
                                             start = 1,
                                             end = nchar(trimmed))
                    } else if (method == 4) {
                        seq[j] <- snpbase
                    }
                }
                
                # If the SNP state is missing, assign NNNNs
                if (is.na(genmat[i, j])) {
                    seq[j] <- "N"
                    if (method == 2) {
                        seq[j] <-
                            stringr::str_pad(
                                seq[j],
                                nchar(trimmed),
                                side = c("right"),
                                pad = "N"
                            )
                    }
                }
            }
            # setTxtProgressBar(pb, i/r)
            
            # Join all the trimmed sequence together into one long 'composite' haplotype
            result <- paste(seq, sep = "", collapse = "")
            # Write the results to file in fastA format
            cat(paste0(">", indNames(x)[i], "_", pop(x)[i], "\n"))
            cat(result, "\n", sep = "")

        }  # Select the next individual and repeat
        
        # Close the output fastA file
        sink()
        
    }
    
    # FLAG SCRIPT END

    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }

    return(invisible(NULL))
    
}
