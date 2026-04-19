#' @name gl.read.vegsnpdb
#' @title Imports vegSNPDB SNP data into dartR and converts it to a genlight object
#' @family io
#'
#' @description
#' Reads SNP genotypes exported from vegSNPDB and converts them into a
#' \code{genlight} object compatible with dartR.
#'
#' The vegSNPDB SNP file is expected to contain locus-level metadata columns
#' (e.g. \code{Chr, Pos, Ref/Alt, ..., MissRate}) followed by one column per
#' individual containing diploid genotypes coded as \code{AA, GG, AG, GA, NN}
#' etc., where codes containing \code{"N"} denote missing data.
#'
#' Genotypes are recoded using the \code{Ref/Alt} column so that:
#' 0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate,
#' NA = missing.
#'
#' Individual IDs are taken from the header row of the genotype columns.
#' In vegSNPDB exports these IDs may contain hyphens/spaces/commas; these are
#' converted to underscores for safe handling in downstream pipelines.
#'
#' If \code{ind.metrics} is provided, it must point to the Excel file
#' (e.g. \code{"Population information.xlsx"}) containing one sheet per crop.
#' The sheet corresponding to \code{vegename} is used to populate individual
#' metadata. Column \code{"Accession ID"} is renamed to \code{id} and
#' \code{"Population"} to \code{pop}. All other column names have spaces
#' replaced by underscores. Individuals in this table are matched to the SNP
#' data via \code{id}.
#'
#' If \code{ind.metrics} is NULL, minimal individual metadata are created with
#' \code{id = indNames(gl)} and \code{pop = "Pop1"}.
#'
#' @param x Name of the vegSNPDB SNP file (csv/tsv/txt, often with .xls
#' extension) [required].
#' @param ind.metrics Path to the Excel file containing per-accession metadata
#' (e.g. \code{"Population information.xlsx"}) [optional but recommended].
#' @param vegename Character. Name of the sheet/table within
#' \code{ind.metrics} to use (e.g. "tomato") [default "tomato"].
#' @param topskip Integer. Number of lines to skip before the header row in
#' the SNP file [default 0].
#' @param lastcolumn Name or index of the last locus-metadata column; genotype
#' columns start immediately after this column [default "MissRate"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, or as set by \code{gl.set.verbosity}].
#'
#' @return A dartR \code{genlight} object containing SNP genotypes and
#' associated locus and individual metadata.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' gl <- gl.read.vegsnpdb(x="Table.xls", ind.metrics="Population information.xlsx",
#'                        vegename="pepper")
#' }
gl.read.vegsnpdb <- function(x,
                             ind.metrics = NULL,
                             vegename    = "pepper",
                             topskip     = 0,
                             lastcolumn  = "MissRate",
                             verbose     = NULL) {
  
  # SET VERBOSITY -------------------------------------------------------------
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START ---------------------------------------------------------
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbose = verbose)
  
  # Helper: make safe labels (hyphens/spaces/commas -> underscores; collapse) --
  .clean_label <- function(z) {
    z <- as.character(z)
    z <- sub("^\ufeff", "", z)                 # strip UTF-8 BOM if present
    z <- gsub("^\\.+", "", z)                  # drop leading "..." (common in exports)
    z <- gsub("[,\\s\\-]+", "_", z)            # commas, whitespace, hyphens -> "_"
    z <- gsub("[^A-Za-z0-9_]", "_", z)         # any other odd chars -> "_"
    z <- gsub("_+", "_", z)                    # collapse multiple underscores
    z <- gsub("^_|_$", "", z)                  # trim leading/trailing underscore
    z
  }
  
  # READ SNP FILE -------------------------------------------------------------
  if (verbose >= 1) {
    cat(report("Reading vegSNPDB SNP file:", x, "\n"))
  }
  
  if (!file.exists(x)) {
    cat(error("Fatal Error: SNP file '", x, "' not found.\n", sep = ""))
    stop()
  }
  
  # topskip is the number of lines to skip BEFORE the header
  skip_lines <- max(as.integer(topskip), 0L)
  
  # Read header line (vegSNPDB export may be mixed-delimiter)
  hdr <- readLines(x, n = skip_lines + 1, encoding = "UTF-8")
  if (length(hdr) < (skip_lines + 1)) {
    cat(error("Fatal Error: File ended before header line could be read.\n"))
    stop()
  }
  hdr <- hdr[skip_lines + 1]
  hdr <- sub("^\ufeff", "", hdr)  # remove BOM if present
  
  # Identify start of comma-separated sample IDs in header
  k <- regexpr(",", hdr)
  
  if (k < 0) {
    # Fallback: standard CSV
    if (verbose >= 2) {
      cat(warn("  Header contains no commas; attempting standard CSV read.\n"))
    }
    df <- utils::read.csv(
      file             = x,
      header           = TRUE,
      stringsAsFactors = FALSE,
      check.names      = FALSE,
      fileEncoding     = "UTF-8-BOM",
      skip             = skip_lines
    )
    
    # Clean all names (safe handling)
    names(df) <- .clean_label(names(df))
    
    # Standardise key locus metadata column names (tolerant to cleaning/case) ----
    .std_locus_names <- function(nm) {
      nm0 <- nm
      nmL <- tolower(nm0)
      
      # helper: rename first match
      .rename_first <- function(from_candidates, to) {
        hit <- which(nmL %in% from_candidates)
        if (length(hit) > 0) {
          nm0[hit[1]] <<- to
          nmL[hit[1]] <<- tolower(to)
        }
      }
      
      # Chr
      .rename_first(c("chr", "chrom", "chromosome"), "Chr")
      
      # Pos (vegSNPDB exports often vary here)
      .rename_first(c("pos", "position", "bp", "basepair", "base_pair"), "Pos")
      
      # Ref/Alt: may have been cleaned to Ref_Alt
      .rename_first(c("ref/alt", "ref_alt", "refalt", "ref_alt_allele"), "Ref/Alt")
      
      nm0
    }
    
    names(df) <- .std_locus_names(names(df))
    
  } else {
    # Parse mixed header into (meta cols) + (sample cols)
    pre  <- substr(hdr, 1, k - 1)               # whitespace-separated meta names
    post <- substr(hdr, k, nchar(hdr))          # comma-separated sample IDs (includes first comma)
    
    pre  <- trimws(gsub("\\s+", " ", pre))
    pre_cols <- strsplit(pre, " ", fixed = TRUE)[[1]]
    
    post_cols <- strsplit(post, ",", fixed = TRUE)[[1]]
    post_cols <- trimws(post_cols)
    post_cols <- post_cols[nzchar(post_cols)]  # drop empties
    
    # Clean labels
    pre_cols  <- .clean_label(pre_cols)
    post_cols <- .clean_label(post_cols)
    
    col_names <- c(pre_cols, post_cols)
    
    # Read the data rows as CSV (skip topskip lines + header line)
    df <- utils::read.csv(
      file             = x,
      header           = FALSE,
      skip             = skip_lines + 1,
      stringsAsFactors = FALSE,
      check.names      = FALSE,
      fileEncoding     = "UTF-8-BOM"
    )
    
    if (ncol(df) != length(col_names)) {
      cat(error(
        "Fatal Error: Parsed ", length(col_names),
        " column names from header but data has ", ncol(df), " columns.\n", sep = ""
      ))
      cat(error("Header preview:\n", hdr, "\n", sep = ""))
      stop()
    }
    
    names(df) <- col_names
  }
  
  if (ncol(df) < 3) {
    cat(error("Fatal Error: vegSNPDB file appears to have too few columns.\n"))
    stop()
  }
  
  # Standardise expected column names after cleaning --------------------------
  # Ref/Alt typically becomes Ref_Alt; map it back so downstream logic works.
  if ("Ref_Alt" %in% names(df) && !("Ref/Alt" %in% names(df))) {
    names(df)[names(df) == "Ref_Alt"] <- "Ref/Alt"
  }
  
  # DETERMINE LAST LOCUS-METADATA COLUMN -------------------------------------
  if (is.character(lastcolumn)) {
    lastcolumn_clean <- .clean_label(lastcolumn)
    # allow for the possibility the file’s column name has been cleaned
    if (!(lastcolumn %in% names(df)) && (lastcolumn_clean %in% names(df))) {
      lastcolumn <- lastcolumn_clean
    }
    if (!(lastcolumn %in% names(df))) {
      cat(error(
        "Fatal Error: lastcolumn '", lastcolumn,
        "' not found in SNP header.\n", sep = ""
      ))
      stop()
    }
    last_meta_col <- which(names(df) == lastcolumn)[1]
  } else if (is.numeric(lastcolumn)) {
    last_meta_col <- as.integer(lastcolumn)
    if (last_meta_col < 1 || last_meta_col >= ncol(df)) {
      cat(error("Fatal Error: lastcolumn index is out of bounds.\n"))
      stop()
    }
  } else {
    cat(error("Fatal Error: lastcolumn must be a column name or numeric index.\n"))
    stop()
  }
  
  if (last_meta_col >= ncol(df)) {
    cat(error(
      "Fatal Error: No genotype columns found after lastcolumn = ",
      lastcolumn, ".\n"
    ))
    stop()
  }
  
  # SPLIT LOCUS METADATA AND GENOTYPE MATRIX ---------------------------------
  locus_meta <- df[, 1:last_meta_col, drop = FALSE]
  geno_char  <- as.matrix(df[, (last_meta_col + 1):ncol(df), drop = FALSE])
  
  n_loci <- nrow(geno_char)
  n_ind  <- ncol(geno_char)
  
  if (verbose >= 1) {
    cat(report("  Detected", n_loci, "loci and", n_ind,
               "individuals in SNP file\n"))
  }
  
  # CHECK REQUIRED LOCUS FIELDS ----------------------------------------------
  required_cols <- c("Chr", "Pos", "Ref/Alt")
  missing_cols  <- setdiff(required_cols, names(locus_meta))
  if (length(missing_cols) > 0) {
    cat(error(
      "Fatal Error: Required columns missing from locus metadata: ",
      paste(missing_cols, collapse = ", "), "\n"
    ))
    stop()
  }
  
  # LOCUS NAMES & BASIC LOCUS METADATA ---------------------------------------
  locus_names <- paste0(locus_meta$Chr, "_", locus_meta$Pos)
  locus_names <- .clean_label(locus_names)  # safe locus labels
  locus_meta$AlleleID <- locus_names
  
  # INDIVIDUAL NAMES FROM HEADER ---------------------------------------------
  ind_names <- .clean_label(colnames(geno_char))
  colnames(geno_char) <- ind_names
  
  # CONVERT GENOTYPES TO 0/1/2/NA --------------------------------------------
  if (verbose >= 1) {
    cat(report("  Converting vegSNPDB genotypes to 0/1/2/NA using Ref/Alt\n"))
  }
  
  refalt_vec <- locus_meta[["Ref/Alt"]]
  if (any(is.na(refalt_vec))) {
    cat(error(
      "Fatal Error: NA values detected in 'Ref/Alt' column; cannot proceed.\n"
    ))
    stop()
  }
  
  refalt_split <- strsplit(refalt_vec, "/", fixed = TRUE)
  bad_ra <- which(lengths(refalt_split) != 2)
  if (length(bad_ra) > 0) {
    cat(error(
      "Fatal Error: 'Ref/Alt' must be formatted like 'A/T'. Bad rows include: ",
      paste(head(bad_ra, 10), collapse = ", "), "\n", sep = ""
    ))
    stop()
  }
  
  ref <- toupper(trimws(vapply(refalt_split, function(z) z[1], character(1))))
  alt <- toupper(trimws(vapply(refalt_split, function(z) z[2], character(1))))
  
  geno_num <- matrix(NA_real_, nrow = n_loci, ncol = n_ind,
                     dimnames = list(locus_names, ind_names))
  
  for (i in seq_len(n_loci)) {
    g <- toupper(trimws(geno_char[i, ]))
    
    # Missing: NA, empty, or containing 'N'
    is_missing <- is.na(g) | g == "" | grepl("N", g, fixed = TRUE)
    g[is_missing] <- NA
    
    r <- ref[i]
    a <- alt[i]
    
    homRef <- paste0(r, r)
    homAlt <- paste0(a, a)
    het1   <- paste0(r, a)
    het2   <- paste0(a, r)
    
    x <- rep(NA_real_, length(g))
    x[g == homRef]           <- 0
    x[g == homAlt]           <- 2
    x[g == het1 | g == het2] <- 1
    
    unknown <- unique(g[
      !(g %in% c(homRef, homAlt, het1, het2)) & !is.na(g)
    ])
    if (length(unknown) > 0 && verbose >= 2) {
      cat(warn(
        "  Warning: Locus ", locus_names[i],
        " has unrecognised genotype codes: ",
        paste(unknown, collapse = ", "),
        " (treated as NA)\n"
      ))
    }
    
    geno_num[i, ] <- x
  }
  
  # TRANSPOSE TO INDIVIDUAL x LOCUS FOR GENLIGHT -----------------------------
  geno_num_t <- t(geno_num)
  
  # CREATE GENLIGHT OBJECT ----------------------------------------------------
  if (verbose >= 1) {
    cat(report("  Creating genlight object\n"))
  }
  
  gl <- new(
    "genlight",
    geno_num_t,
    ploidy    = 2,
    loc.names = locus_names,
    ind.names = ind_names
  )
  
  # STORE CHROMOSOME AND POSITION (SLOTS) ------------------------------------
  gl@chromosome <- factor(locus_meta$Chr)
  
  gl@position <- suppressWarnings(as.integer(locus_meta$Pos))
  if (anyNA(gl@position)) {
    bad <- which(is.na(gl@position))[1]
    stop(
      "Fatal Error: Non-numeric/NA positions detected in 'Pos' column; e.g. locus ",
      locus_names[bad]
    )
  }
  
  # LOCUS METRICS -------------------------------------------------------------
  gl@other$loc.metrics <- locus_meta
  
  # INDIVIDUAL METRICS FROM 'Population information...' ----------------------
  if (!is.null(ind.metrics)) {
    
    if (!file.exists(ind.metrics)) {
      cat(error(
        "Fatal Error: Individual metadata file '", ind.metrics,
        "' not found.\n", sep = ""
      ))
      stop()
    }
    
    if (!requireNamespace("readxl", quietly = TRUE)) {
      cat(error(
        "Fatal Error: Package 'readxl' is required to read Excel metadata.\n"
      ))
      stop()
    }
    
    sheets <- readxl::excel_sheets(ind.metrics)
    if (!(vegename %in% sheets)) {
      cat(error(
        "Fatal Error: vegename '", vegename,
        "' not found in '", basename(ind.metrics),
        "'. vegename must be one of ",
        paste(sheets, collapse = ", "),
        ".\n", sep = ""
      ))
      stop()
    }
    
    if (verbose >= 1) {
      cat(report(
        "  Reading individual metadata from sheet '", vegename,
        "' in ", ind.metrics, "\n", sep = ""
      ))
    }
    
    # Your file has a one-line title row, then the header row.
    meta_raw <- readxl::read_excel(
      ind.metrics,
      sheet    = vegename,
      skip     = 1
    )
    
    # Standardise column names: spaces/hyphens/commas -> underscores
    cn <- colnames(meta_raw)
    cn <- .clean_label(cn)
    
    # Rename specific columns (after cleaning)
    cn[cn == "Accession_ID"] <- "id"
    cn[cn == "AccessionID"]  <- "id"
    cn[cn == "Population"]   <- "pop"
    
    colnames(meta_raw) <- cn
    
    if (!("id" %in% cn)) {
      cat(error(
        "Fatal Error: Individual metadata sheet '", vegename,
        "' must contain column 'Accession ID' (mapped to 'id').\n", sep = ""
      ))
      stop()
    }
    if (!("pop" %in% cn)) {
      cat(error(
        "Fatal Error: Individual metadata sheet '", vegename,
        "' must contain column 'Population' (mapped to 'pop').\n", sep = ""
      ))
      stop()
    }
    
    meta_df <- as.data.frame(meta_raw, stringsAsFactors = FALSE)
    
    # Clean the id values to match the SNP individual IDs if needed
    meta_df$id <- .clean_label(meta_df$id)
    
    # MATCH SNP INDIVIDUALS TO METADATA VIA id -------------------------------
    mm <- match(ind_names, meta_df$id)
    if (any(is.na(mm))) {
      missing_ids <- ind_names[is.na(mm)]
      cat(error(
        "Fatal Error: The following individuals from SNP data were not found ",
        "in metadata (id column) for vegename '", vegename, "': ",
        paste(missing_ids, collapse = ", "), "\n", sep = ""
      ))
      stop()
    }
    
    meta_df <- meta_df[mm, , drop = FALSE]
    
    gl@other$ind.metrics <- meta_df
    pop(gl) <- factor(meta_df$pop)
    
  } else {
    # DEFAULT: POP1 FOR ALL INDIVIDUALS --------------------------------------
    if (verbose >= 1) {
      cat(report("  No individual metadata file provided; setting pop = 'Pop1'\n"))
    }
    ind_meta <- data.frame(
      id  = ind_names,
      pop = factor("Pop1", levels = "Pop1"),
      stringsAsFactors = FALSE
    )
    gl@other$ind.metrics <- ind_meta
    pop(gl) <- ind_meta$pop
  }
  
  # COMPLIANCE CHECK ----------------------------------------------------------
  gl <- gl.compliance.check(gl, verbose = verbose)
  
  # HISTORY -------------------------------------------------------------------
  gl@other$history <- list()
  gl@other$history[[1]] <- match.call()
  
  # SUMMARY -------------------------------------------------------------------
  if (verbose >= 3) {
    cat("\nSummary of vegSNPDB SNP dataset\n")
    cat("  No. of loci:       ", nLoc(gl), "\n")
    cat("  No. of individuals:", nInd(gl), "\n")
    cat("  No. of populations:", nPop(gl), "\n\n")
  }
  
  # FLAG SCRIPT END -----------------------------------------------------------
  if (verbose > 0) {
    cat(report(paste("Completed:", funname, "\n")))
  }
  
  return(gl)
}
