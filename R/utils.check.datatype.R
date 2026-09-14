#' @name utils.check.datatype
#' @title Utility function to check the class of an object passed to a function
#' @family utilities
#' 
#' @description 
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.

#' @param x Name of the genlight object, dist matrix, data matrix, glPCA, or
#' fixed difference list (fd) [required].
#' @param accept Vector containing the classes of objects that are to be
#' accepted. Matching is exact and case-sensitive. The entries 'genlight'
#' and 'dartR' accept both genotype datatypes ('SNP' and 'SilicoDArT')
#' [default c('genlight','SNP','SilicoDArT','dartR')].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#' 
#' @details
#' Most functions require access to a genlight object, dist matrix, data matrix
#'  or fixed difference list (fd), and this function checks that a genlight
#'  object or one of the above has been passed, whether the genlight object is a
#'   SNP dataset or a SilicoDArT object, and reports back if verbosity is >=2.
#'   
#' This function checks the class of passed object and sets the datatype to
#' 'SNP', 'SilicoDArT', 'dist', 'matrix', 'glPca', 'fd', 'data.frame',
#' 'list', or class(x)[1] as appropriate.

#' A genlight object is classified from its ploidy slot alone (uniformly 1,
#' 'SilicoDArT'; uniformly 2, 'SNP'; mixed or other ploidy is treated as
#' 'SNP', with a warning) -- genotype values, dartR flags and locus metrics
#' are not consulted, with one exception: an object of uniform ploidy 2
#' whose non-missing genotypes are all 0 or 1, and which carries no SNP
#' metadata (an empty loc.all slot and no SNP or SnpPosition locus
#' metrics), is rejected as presence/absence (SilicoDArT) content
#' mislabelled as SNP, since it would otherwise pass into dosage-based
#' (0/1/2) arithmetic. A clean SNP subset that merely lacks the
#' homozygous-alternate class is not affected: its loc.all slot vouches
#' for the label. Content-level validation beyond that check is the job
#' of \code{gl.compliance.check}.

#' Note also that this function checks to see if there are individuals or loci
#' scored as all missing (NA) and if so, issues the user with a warning
#' (verbose >= 2).

#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' datatype <- utils.check.datatype(testset.gl)
#' datatype <- utils.check.datatype(as.matrix(testset.gl),accept='matrix')
#' fd <- gl.fixed.diff(testset.gl)
#' datatype <- utils.check.datatype(fd,accept='fd')
#'  
#' @return datatype (returned invisibly): 'SNP' for SNP data, 'SilicoDArT' for
#' P/A data, 'dist' for a distance matrix, 'matrix' for a data matrix, 'glPca'
#' for an ordination file, 'fd' for a fixed difference object, 'data.frame'
#' for a data frame, 'list' for any other list, or class(x)[1].
#' @export

utils.check.datatype <- function(x,
                                 accept = c("genlight", "SNP", "SilicoDArT","dartR"),
                                 verbose = NULL) {
    #### SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    if (is(x, "dartR") || is(x, "genlight")) {
        if (is.null(ploidy(x))) {
            stop(
                error(
                    "Fatal Error: ploidy not set in the genlight object, run gl <- gl.compliance.check(gl)\n"
                )
            )
        }
        if (verbose >= 2) {
            cat(report("  Processing genlight object"))
        }
        if (all(ploidy(x) == 1)) {
            if (verbose >= 2) {
                cat(report(" with Presence/Absence (SilicoDArT) data\n"))
            }
            datatype <- "SilicoDArT"
        } else if (all(ploidy(x) == 2)) {
            if (verbose >= 2) {
                cat(report(" with SNP data\n"))
            }
            datatype <- "SNP"
            # Content-vs-ploidy consistency: an object whose non-missing
            # genotypes are all 0 or 1 is presence/absence (SilicoDArT)
            # content mislabelled with ploidy 2 -- left unchecked it
            # passes accept = 'SNP' gates into dosage-based (0/1/2)
            # arithmetic. The scan exits at the first genotype of 2, so
            # clean SNP data pay for a single individual (or a single
            # column block when FBM-backed); an object with no non-NA
            # genotypes at all is left to the all-NA warnings below.
            # (a zero-locus object has no content to contradict; note
            # as.integer() on a zero-locus SNPbin returns a spurious 0)
            seen.value <- FALSE
            seen.dosage <- FALSE
            if (nLoc(x) == 0) {
                # no genotypes to scan
            } else if (methods::.hasSlot(x, "fbm") && !is.null(x@fbm)) {
                p <- ncol(x@fbm)
                start <- 1L
                while (start <= p && !seen.dosage) {
                    idx <- start:min(start + 1023L, p)
                    blk <- x@fbm[, idx, drop = FALSE]
                    if (!seen.value && any(!is.na(blk))) {
                        seen.value <- TRUE
                    }
                    if (any(blk > 1, na.rm = TRUE)) {
                        seen.dosage <- TRUE
                    }
                    start <- start + 1024L
                }
            } else {
                for (ind in x@gen) {
                    v <- as.integer(ind)
                    if (!seen.value && any(!is.na(v))) {
                        seen.value <- TRUE
                    }
                    if (any(v > 1, na.rm = TRUE)) {
                        seen.dosage <- TRUE
                        break
                    }
                }
            }
            if (seen.value && !seen.dosage) {
                # A clean SNP subset can lack the homozygous-alternate
                # class entirely (a protected use case -- see
                # test-gl.report.basics.R, test-gl.smearplot.R), so
                # 0/1-only content alone is not proof of mislabelling.
                # SNP-level metadata vouches for the label: the loc.all
                # slot (the two alleles per locus -- definitionally SNP,
                # guaranteed by gl.compliance.check for SNP data) or the
                # SNP/SnpPosition locus metrics from a DArT SNP read.
                # SilicoDArT objects carry neither (verified: testset.gs
                # has an empty loc.all and lacks both columns even after
                # SNP-style metric recalculation).
                snp.evidence <- length(x@loc.all) > 0
                if (!snp.evidence && !is.null(x@other$loc.metrics)) {
                    snp.evidence <- any(c("SNP", "SnpPosition") %in%
                                            names(x@other$loc.metrics))
                }
                if (!snp.evidence) {
                    stop(
                        error(
                            "Fatal Error: object has ploidy 2 (SNP) but all non-missing genotypes are scored 0 or 1 and no SNP metadata (loc.all, SNP or SnpPosition locus metrics) is present; ploidy slot and genotype content disagree. If the data are Tag P/A (SilicoDArT), reassign ploidy with ploidy(gl) <- rep(1, nInd(gl)) or re-read the data\n"
                        )
                    )
                }
            }
        } else {
          # Mixed or non-diploid ploidy (e.g. polyploid data) is treated
          # as SNP data, with notice -- the notice affects how results
          # are computed downstream, so it prints from verbose 1
          if (verbose >= 2) {
            cat(report(" with SNP data\n"))
          }
          if (verbose >= 1) {
            cat(warn(
              "  Warning: ploidy is not uniformly 2; treating as SNP data\n"
            ))
          }
          datatype <- "SNP"
        }
        # Check for individuals or loci scoring all missing values (NA)
        # (a direct single-pass check; running a full gl.filter.allna
        # here cost ~50 ms on every function entry)
        if (verbose > 1) {
            if (methods::.hasSlot(x, "fbm") && !is.null(x@fbm)) {
                # FBM-backed: scan in column blocks via the FBM accessor
                # rather than materialising the full matrix (DAT6)
                n <- nrow(x@fbm)
                p <- ncol(x@fbm)
                loc.allna <- FALSE
                ind.scored <- integer(n)
                start <- 1L
                while (start <= p) {
                    idx <- start:min(start + 1023L, p)
                    blk <- !is.na(x@fbm[, idx, drop = FALSE])
                    if (!is.matrix(blk)) {
                        blk <- matrix(blk, nrow = n)
                    }
                    if (any(colSums(blk) == 0)) {
                        loc.allna <- TRUE
                    }
                    ind.scored <- ind.scored + rowSums(blk)
                    start <- start + 1024L
                }
                ind.allna <- any(ind.scored == 0)
            } else {
                namat <- is.na(as.matrix(x))
                loc.allna <- any(colSums(!namat) == 0)
                ind.allna <- any(rowSums(!namat) == 0)
            }
            if (loc.allna) {
                cat(
                    warn(
                        "  Warning: data include loci that are scored NA across all individuals.\n  Consider filtering using gl <- gl.filter.allna(gl)\n"
                    )
                )
            }
            if (ind.allna) {
                cat(
                    warn(
                        "  Warning: data include individuals that are scored NA across all loci.\n  Consider filtering using gl <- gl.filter.allna(gl)\n"
                    )
                )
            }
        }
    } else if (is(x, "fd")) {
        if (is(x$gl, "genlight")) {
            # if(is.null(ploidy(x$gl))){ stop(error('Fatal Error: ploidy not set in the genlight object, run gl <-
            # gl.compliance.check(gl)\n')) }
            if (verbose >= 2) {
                cat(report("  Processing a fixed difference (fd) object"))
            }
            if (verbose >= 2) {
                if (all(ploidy(x$gl) == 1)) {
                    cat(report(
                        " with Presence/Absence (SilicoDArT) data\n"
                    ))
                } else {
                    cat(report(" with SNP data\n"))
                }
            }
        } else {
            stop(
                error(
                    "Fatal Error: Fixed Difference object expected! Check format of object\n"
                )
            )
        }
        datatype <- "fd"
    } else if (is(x, "dist")) {
        if (verbose >= 2) {
            cat(report("  Processing a distance matrix\n"))
        }
        datatype <- "dist"
    } else if (is(x, "matrix")) {
        if (verbose >= 2) {
            cat(report("  Processing a data matrix\n"))
        }
        datatype <- "matrix"
    } else if (is(x, "glPca")) {
        if (verbose >= 2) {
            cat(report("  Processing an ordination file (glPca)\n"))
        }
        datatype <- "glPca"
    } else if (is.data.frame(x)) {
        if (verbose >= 2) {
            cat(report("  Processing a data frame\n"))
        }
        datatype <- "data.frame"
    } else if (is.list(x)) {
        if (verbose >= 2) {
            cat(report("  Processing a list\n"))
        }
        datatype <- "list"
    } else {
        if (verbose >= 1) {
            cat(warn("  Warning: Found object of class", class(x)[1], "\n"))
        }
        datatype <- class(x)[1]
    }

    #### CHECK WHETHER TO THROW AN ERROR ####

    # 'genlight' or 'dartR' in accept admits both genotype datatypes,
    # unless a specific genotype datatype is also listed (in which case
    # the specific listing governs, e.g. c('genlight','SNP') is SNP-only)
    if (datatype %in% c("SNP", "SilicoDArT") &&
        any(c("genlight", "dartR") %in% accept) &&
        !any(c("SNP", "SilicoDArT") %in% accept)) {
        accept <- union(accept, c("SNP", "SilicoDArT"))
    }

    if (!(datatype %in% accept)) {
        stop(
            error(
                "Fatal Error: inappropriate object passed to function, found",
                datatype,
                "expecting",
                paste0(paste(accept, collapse = " or "), "\n")
            )
        )
    }
    
    invisible(datatype)
}

