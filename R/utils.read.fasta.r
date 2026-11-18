#' @name utils.read.fasta
#' @title An internal script to read a fastA file into a genlight object
#' @family utilities
#' 
#' @description 
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.
#'
#' @param file Name of the fastA file [required]
#' @param parallel Switch to deactivate parallel version. It might not be worth
#' to run it parallel most of the times [default FALSE]
#' @param n.cores Number of cores to use if parallel = TRUE; if NULL, all
#' the available cores are used [default NULL]
#' @param verbose Verbosity: 0, silent, fatal errors only; 1, flag function
#' begin and end; 2, progress log; 3, progress and results summary; 5, full
#' report [default 2 or as specified using gl.set.verbosity].
#'  
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})

# @export
#' @return The resultant genlight object

utils.read.fasta <-  function(file,
                              parallel = FALSE,
                              n.cores = NULL,
                              verbose = NULL) {

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  if(verbose >= 2){
  cat(report("  Reading",basename(file),"\n"))
  }

  # Resolve the parallel settings up front: n.cores = NULL means all the
  # available cores, and fork-based mclapply is unavailable on Windows,
  # where the read falls back to serial
  if (parallel) {
    if (is.null(n.cores)) {
      n.cores <- parallel::detectCores()
    }
    if (.Platform$OS.type == "windows" && n.cores > 1) {
      if (verbose >= 2) {
        cat(warn(
          "  Warning: fork-based parallelism is not available on Windows; reading serially\n"
        ))
      }
      parallel <- FALSE
    }
  }

  txt <- scan(file,
              what = "character",
              sep = "\n",
              quiet = TRUE)

  # Normalise softmasked (lowercase) bases once on read; a lowercase base
  # must not register as a distinct allele. Headers keep their case
  is_header <- grepl("^>", txt)
  txt[!is_header] <- toupper(txt[!is_header])

  nb.ind <- sum(is_header)
  # find individuals' labels
  IND.LAB <- sub(">", "", txt[is_header])
  # split per individuals (two lines per record: header + sequence)
  txt <- split(txt, rep(1:nb.ind, each = 2))
  if (parallel) {
    # each genome -> one vector
    seqs <-
      parallel::mclapply(txt, function(e)
        strsplit(paste(e[-1], collapse = ""), split = "")[[1]],
        mc.cores = n.cores, mc.silent = TRUE, mc.cleanup =
          TRUE, mc.preschedule = FALSE)
  } else {
    # each genome -> one vector
    seqs <-
      lapply(txt, function(e)
        strsplit(paste(e[-1], collapse = ""), split = "")[[1]])
  }

  # All sequences must be aligned to the same length; recycling a short
  # record would silently misalign every subsequent locus
  NLOC <- length(seqs[[1]])
  seq.lengths <- lengths(seqs)
  if (any(seq.lengths != NLOC)) {
    stop(error(
      "Fatal Error: sequences in", basename(file),
      "do not all have the same length; check record(s):",
      paste(unique(IND.LAB[seq.lengths != NLOC]), collapse = ", "), "\n"
    ))
  }

  seq_mat <- matrix(unlist(seqs), byrow = TRUE, nrow = nb.ind)

  # Genotype expansion table: homozygotes and the one-letter IUPAC
  # heterozygotes; anything else (N, -, V, H, D, B, unrecognized) is
  # missing data and contributes no alleles
  geno_expand <- list(
    A = c("A", "A"), C = c("C", "C"), G = c("G", "G"), T = c("T", "T"),
    M = c("A", "C"), R = c("A", "G"), W = c("A", "T"),
    S = c("C", "G"), Y = c("C", "T"), K = c("G", "T")
  )

  ## POOL contains the alleles of each position, derived from the expanded
  ## genotype codes so that a heterozygote contributes both its alleles
  POOL <- lapply(seq_len(ncol(seq_mat)), function(j) {
    sort(unique(unlist(geno_expand[unique(seq_mat[, j])])))
  })

  nb.alleles <- lengths(POOL)
  snp.posi <- which(nb.alleles == 2)
  multi.posi <- which(nb.alleles > 2)

  # Columns whose allele pool exceeds two letters cannot be scored as
  # biallelic SNPs and are skipped; the message affects what the output
  # contains, so it prints from verbose = 1 up
  if (length(multi.posi) > 0 && verbose >= 1) {
    cat(important(
      "  SNP positions", paste(multi.posi, collapse = " "),
      "from file", basename(file),
      "have more than 2 alleles. They are skipped\n"
    ))
  }

  if (length(snp.posi) == 0) {
    if (verbose >= 1) {
      cat(warn("  No polymorphism in the alignment", basename(file), "\n"))
    }
    return(NULL)
  }

  # Classify genotypes per retained column. The reference allele is the
  # allele (not genotype class) with the highest count over the expanded
  # genotypes, ties resolving alphabetically; the dosage is the count of
  # the alternate allele (0 hom ref, 1 het, 2 hom alt), NA for missing
  snp_calls <- lapply(snp.posi, function(j) {
    expanded <- geno_expand[seq_mat[, j]]
    counts <- table(factor(unlist(expanded), levels = POOL[[j]]))
    ref <- names(counts)[which.max(counts)]
    alt <- setdiff(POOL[[j]], ref)
    dosage <- vapply(seq_len(nb.ind), function(i) {
      g <- expanded[[i]]
      if (is.null(g)) NA_integer_ else sum(g == alt)
    }, integer(1))
    list(dosage = dosage, loc.all = paste0(ref, "/", alt))
  })

  snp_mat <- do.call(cbind, lapply(snp_calls, "[[", "dosage"))
  loc.all <- vapply(snp_calls, "[[", character(1), "loc.all")

  res <- apply(snp_mat, 1, function(e)
      new(
        "SNPbin", snp = e, ploidy = 2L
      ))

  res <- new("genlight", res, ploidy = 2L)

  indNames(res) <- IND.LAB
  alleles(res) <- loc.all
  locNames(res) <-  paste0(sub("\\..*", "", basename(file)), "_", snp.posi)

  return(res)

}

merge_gl_fasta <- function(gl_list, 
                           parallel = FALSE,
                           verbose = verbose) {
  
  if(verbose >= 2){
    cat(report("  Merging files...\n"))
  }
  
 mono_file <- unlist(lapply(gl_list,function(x){class(x)[[1]]}))
 
 mono_file <- which(mono_file=="NULL")
 
 if(length(mono_file)>0){
   gl_list <- gl_list[-mono_file]
 }
  
  if(length(gl_list)==1){
    res <- gl_list[[1]]
  }else{

  # Duplicate individual names within a file would be joined many-to-one
  # by merge(), silently copying one record's genotypes onto several
  # rows; refuse them before any join happens
  for (y in gl_list) {
    dups <- unique(indNames(y)[duplicated(indNames(y))])
    if (length(dups) > 0) {
      stop(error(
        "Fatal Error: duplicate individual name(s)",
        paste(dups, collapse = ", "),
        "in file", sub("_[0-9]+$", "", locNames(y)[1]),
        "-- rename the records before merging fasta files\n"
      ))
    }
  }

  matrix_temp <- lapply(gl_list, function(y) {
    return(as.data.frame(cbind(ind_names = indNames(y), as.matrix(y))))
  })
  
  loc_names <- lapply(matrix_temp,function(x){
  colnames(x)[-1]  
  })

  loc_names <- Reduce("c",loc_names)
  
  loc_all <- lapply(gl_list, function(y) {
    return(y$loc.all)
  })
  
  loc_all <- Reduce("c",loc_all)
  
  gl_temp <-
    Reduce(function(x, y) {
      merge(x, y, by = "ind_names", all = TRUE)
    }, matrix_temp)
  
  res_temp <-  matrix2gen(gl_temp[, 2:ncol(gl_temp)], parallel = parallel)
  
  res <- new("genlight", res_temp, ploidy = 2L)
  
  res$loc.names <- loc_names
  alleles(res) <- loc_all
  res$ind.names <- gl_temp$ind_names
  }
  
  return(res)
  
}
