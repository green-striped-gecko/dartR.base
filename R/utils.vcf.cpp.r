## Internal helpers of gl.read.vcf
##
## .vcf_fix_info() builds the locus information of gl.read.vcf from the
## fixed columns of a vcf (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, as
## vcfR's @fix: a character matrix with "." read as NA). It is shared by the
## vcfR and cpp parsers, so both give the same QUAL/FILTER/INFO metrics.
## Records with more than one ALT allele are dropped, as gl.read.vcf always
## did. Returns the info table, alleles, chromosome and position of the
## kept records and their row indices (keep).
.vcf_fix_info <- function(fix, pos, loc.all) {
  chrom <- fix[, "CHROM"]
  info_tmp_1 <- fix[, 6:7]
  info_tmp_2 <- fix[, "INFO"]
  if (any(is.na(info_tmp_2[1]) | is.na(info_tmp_1[1])) == TRUE) {
    info <- info_tmp_1
    colnames(info) <- c("QUAL", "FILTER")
  } else {
    # Parse each record's INFO string into key=value pairs and fill by
    # key, so records whose keys differ in order or presence still land
    # in the right columns (F2). Flag-type keys (no '=') are recorded as
    # "TRUE"; keys absent from a record are NA.
    info_list <- lapply(stringr::str_split(info_tmp_2, ";"), function(rec) {
      rec <- rec[rec != "" & rec != "."]
      has_eq <- grepl("=", rec, fixed = TRUE)
      keys <- ifelse(has_eq, sub("=.*$", "", rec), rec)
      vals <- ifelse(has_eq, sub("^[^=]*=", "", rec), "TRUE")
      stats::setNames(vals, keys)
    })
    keys_all <- unique(unlist(lapply(info_list, names)))
    if (length(keys_all) == 0) {
      info <- info_tmp_1
      colnames(info) <- c("QUAL", "FILTER")
    } else {
      info_tmp_2 <- as.data.frame(
        do.call(rbind, lapply(info_list, function(rec) unname(rec[keys_all]))),
        stringsAsFactors = FALSE)
      info <- cbind(info_tmp_1, info_tmp_2)
      colnames(info) <- c("QUAL", "FILTER", keys_all)
    }
  }

  keep <- seq_len(nrow(fix))
  # identify which SNPs have more than 2 alleles
  if ("AC" %in% colnames(info)) {
    more_alleles <- grep(",", info$AC)
    if (length(more_alleles) != 0) {
      info <- info[-more_alleles, ]
      # Numerify only the columns that parse as numeric; coercing every
      # column destroyed FILTER ("PASS" -> NA) and any character INFO
      # field whenever a multi-allelic record was present (F4).
      numerifiable <- vapply(info, function(col) {
        col <- as.character(col)
        all(is.na(col) | !is.na(suppressWarnings(as.numeric(col))))
      }, logical(1))
      info[numerifiable] <- lapply(info[numerifiable],
                                   function(col) as.numeric(as.character(col)))
      keep <- keep[-more_alleles]
    }
  } else {
    more_alleles <- grep(pattern = ",", fix[, "ALT"])
    if (length(more_alleles) > 0) {
      info <- info[-more_alleles, ]
      keep <- keep[-more_alleles]
    }
  }
  list(info = info, keep = keep, loc.all = loc.all[keep],
       chromosome = as.factor(chrom[keep]), position = pos[keep])
}

## Fixed columns of a vcf (CHROM to INFO), read in chunks of lines so the
## genotypes are never held in memory: a character matrix with "." as NA,
## as vcfR's @fix. Also reports whether any genotype call (GT, the first
## FORMAT field) is haploid (one allele, e.g. "1") or polyploid (more than
## two alleles), which the cpp parser cannot code as gl.read.vcf does.
.vcf_read_fixed <- function(vcffile, chunk.bytes = 32e6) {
  con <- file(vcffile, open = "r")
  on.exit(close(con), add = TRUE)
  fix <- list()
  n.haploid <- 0
  n.polyploid <- 0
  k <- 1000L
  repeat {
    lines <- readLines(con, n = k, warn = FALSE)
    if (length(lines) == 0) break
    lines <- lines[!startsWith(lines, "#")]
    if (length(lines) == 0) next
    # adapt the number of lines per chunk to the line length
    k <- max(1L, as.integer(chunk.bytes / max(1, mean(nchar(lines[1:min(10, length(lines))])))))
    head8 <- sub("^(([^\t]*\t){7}[^\t]*).*$", "\\1", lines)
    f <- do.call(rbind, strsplit(head8, "\t", fixed = TRUE))
    fix[[length(fix) + 1L]] <- f
    # genotype fields: everything after the 9th column; GT is the first
    # sub-field of each sample
    smp <- sub("^([^\t]*\t){9}", "", lines)
    gt <- gsub("(^|\t)([^:\t]*)[^\t]*", "\\1\\2", smp)
    calls <- unlist(strsplit(gt, "\t", fixed = TRUE), use.names = FALSE)
    calls <- calls[calls != "." & calls != "./." & calls != ".|."]
    nsep <- nchar(gsub("[^/|]", "", calls))
    n.haploid <- n.haploid + sum(nsep == 0)
    n.polyploid <- n.polyploid + sum(nsep > 1)
  }
  fix <- do.call(rbind, fix)
  if (is.null(fix)) {
    fix <- matrix(character(0), 0, 8)
  }
  colnames(fix) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                     "INFO")
  fix[fix == "."] <- NA
  list(fix = fix, n.haploid = n.haploid, n.polyploid = n.polyploid)
}

## The vcf read with tidypopgen's C++ parser (gl.read.vcf(parser = "cpp")).
## Genotypes (count of the ALT allele) are read by tidypopgen into its own
## file-backed matrix and copied in blocks of loci into a dartR FBM; the
## fixed columns are read in chunks (.vcf_read_fixed). The records
## tidypopgen kept (single-base REF and ALT) are matched to the fixed
## columns by CHROM, POS, REF and ALT in file order, so locus names (vcfR's
## addID: ID, or CHROM_POS when ID is missing), alleles, chromosome,
## position and the QUAL/FILTER/INFO metrics are those of the vcfR route.
## Returns the dartR object (file-backed) and its info table.
.vcf_read_cpp <- function(vcffile, verbose = 2) {
  fx <- .vcf_read_fixed(vcffile)
  if (fx$n.haploid > 0 || fx$n.polyploid > 0) {
    stop(error("  Fatal Error: the vcf has", fx$n.haploid, "haploid and",
               fx$n.polyploid, "polyploid genotype calls. parser = \"cpp\"",
               "reads diploid calls only (tidypopgen codes a haploid call",
               "like a heterozygote); use parser = \"vcfR\"\n"))
  }
  fix <- fx$fix
  pos <- as.integer(fix[, "POS"])
  loc.all <- paste0(fix[, "REF"], "/", fix[, "ALT"])
  fi <- .vcf_fix_info(fix, pos = pos, loc.all = loc.all)

  if (verbose >= 2) {
    cat(report("  Reading the genotypes with tidypopgen's C++ parser\n"))
  }
  bk <- utils.fbm.tempfile("tpg_")
  on.exit(.fbm_try_delete(paste0(bk, ".bk")), add = TRUE)
  on.exit(unlink(paste0(bk, ".rds")), add = TRUE)
  gt <- tidypopgen::gen_tibble(vcffile, backingfile = bk, parser = "cpp",
                               quiet = TRUE)
  loci <- tidypopgen::show_loci(gt)
  ids <- as.character(gt$id)
  n <- length(ids)
  p <- nrow(loci)

  # the kept records, matched in file order (n-th occurrence of a key)
  occ <- function(k) paste(k, stats::ave(seq_along(k), k, FUN = seq_along))
  key.fix <- occ(paste(fix[, "CHROM"], pos, fix[, "REF"], fix[, "ALT"]))
  key.tpg <- occ(paste(loci$chromosome, loci$position, loci$allele_ref,
                       loci$allele_alt))
  rows <- match(key.tpg, key.fix)
  if (anyNA(rows) || is.unsorted(rows, strictly = TRUE) ||
      anyNA(match(rows, fi$keep))) {
    stop(error("  Fatal Error: the loci read by tidypopgen could not be",
               "matched to the records of the vcf; use parser = \"vcfR\"\n"))
  }
  n.skip <- length(fi$keep) - p
  if (n.skip > 0 && verbose >= 1) {
    cat(warn("  Warning:", n.skip, "biallelic records that are not",
             "single-base SNPs (e.g. indels) were skipped by parser =",
             "\"cpp\"; use parser = \"vcfR\" to keep them\n"))
  }

  # genotypes in blocks of loci into a dartR FBM (missing as code 3)
  G <- utils.fbm.new(n, p, prefix = "vcf_")
  for (cols in .fbm_col_blocks(p, n)) {
    g <- tidypopgen::show_genotypes(gt, loci_indices = cols)
    g[is.na(g)] <- 3
    G[, cols] <- g
  }
  rm(gt)

  id <- fix[rows, "ID"]
  id[is.na(id)] <- paste(fix[rows, "CHROM"], pos[rows], sep = "_")[is.na(id)]
  k <- match(rows, fi$keep)
  info <- fi$info[k, , drop = FALSE]
  x <- new("dartR",
           fbm = G,
           ploidy = rep(2L, n),
           ind.names = ids,
           loc.names = id,
           loc.all = fi$loc.all[k],
           chromosome = droplevels(fi$chromosome[k]),
           position = fi$position[k])
  list(x = x, info = info)
}
