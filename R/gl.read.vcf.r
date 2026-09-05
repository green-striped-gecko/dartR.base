#' Converts a vcf file into a genlight object
#'
#' This function needs package vcfR, please install it. 
#' @param vcffile A vcf file (works only for diploid data) [required].
#' @param ind.metafile Optional file in csv format with metadata for each
#' individual (see details for explanation) [default NULL].
#' @param mode "genotype" all heterozygous sites will be coded as 1 regardless ploidy level, 
#' dosage: sites will be codes as copy number of alternate allele [default genotype]
#' @param fbm Logical, if TRUE the dartR object will contain a file-backed matrix. 
#' this is important for large data sets that do not fit into RAM.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' The ind.metadata file needs to have very specific headings. First a heading
#' called id. Here the ids have to match the ids in the dartR object.
#' The following column headings are optional.
#' pop: specifies the population membership of each individual. lat and lon
#' specify spatial coordinates (in decimal degrees WGS1984 format). Additional
#' columns with individual metadata can be imported (e.g. age, gender).
#' Individuals in the SNP data that are absent from the ind.metafile are
#' retained in the returned object with NA metadata; a warning lists them
#' at verbose >= 1.
#'
#' Genotypes are coded as the count of the ALT allele: 0 = homozygous
#' REF, 2 (or the copy number in dosage mode) = homozygous ALT. Note that
#' this is the opposite orientation to \code{\link{gl.read.PLINK}}, which
#' counts allele.2 (the PLINK 1.x major allele).
#'
#' In "genotype" mode the returned object is coded as diploid; if the vcf
#' contains haploid or polyploid genotype calls a warning is issued at
#' verbose >= 1 (haploid calls are recoded as diploid homozygotes,
#' polyploid heterozygous calls as 1). In "dosage" mode the ploidy of the
#' object is set to the data's maximum copy number of the alternate
#' allele, uniform across individuals.
#' Note also that this function checks to see if there are input of mode, missing input of mode
#' will issue the user with an error.
#' Please carefully check the data if "dosage" mode is used.
#' @return A genlight object.
#' @export
#' @author Bernd Gruber, Ching Ching Lau (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#' # read in vcf and convert to format as DArT data
#' obj <- gl.read.vcf(system.file('extdata/test.vcf',package='dartR'),
#' ind.metafile="metafile.csv")
#' # read in vcf and convert to format as dosage
#' obj <- gl.read.vcf(system.file('extdata/test.vcf',package='dartR'),
#' ind.metafile="metafile.csv",mode="dosage")
#' }

gl.read.vcf <- function(vcffile,
                        ind.metafile = NULL,
                        mode="genotype",
                        fbm=FALSE,
                        verbose = NULL) {

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jackson",
                   verbose = verbose)
  
  x <- NULL
  
  pkg <- "vcfR"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  } 
  
  vcf <- vcfR::read.vcfR(file = vcffile, verbose = verbose)
  myRef <- vcfR::getREF(vcf)
  myAlt <- vcfR::getALT(vcf)
  chrom <- vcfR::getCHROM(vcf)
  pos <- vcfR::getPOS(vcf)
  loc.all <- paste0(myRef,"/",myAlt)

  # Enforce the documented diploid contract for genotype mode (F6): warn
  # when the GT calls are not diploid -- haploid calls are recoded as
  # diploid homozygotes, polyploid heterozygous calls as 1, and the
  # object is coded diploid.
  if (mode == "genotype" & verbose >= 1) {
    gt_arity <- vcfR::extract.gt(vcf)
    gt_arity <- gt_arity[!is.na(gt_arity)]
    if (length(gt_arity) > 0) {
      gt_arity <- nchar(gsub("[^/|]", "", gt_arity)) + 1L
      if (any(gt_arity != 2)) {
        cat(warn(
          "  Warning: non-diploid genotype calls detected. In 'genotype' mode haploid calls are recoded as diploid homozygotes and polyploid heterozygous calls are coded 1; the returned object is coded as diploid. Use mode='dosage' to retain allele copy number.\n"
        ))
      }
    }
  }

  x <- utils.vcfr2genlight.polyploid(x=vcf, mode2=mode)
  
  # adding SNP information from VCF
  info_tmp_1 <- vcf@fix[,6:7]
  info_tmp_2 <- vcfR::getINFO(vcf)
  if(any(is.na(info_tmp_2[1]) | is.na(info_tmp_1[1]))==TRUE){
    info <- info_tmp_1
    colnames(info) <- c("QUAL","FILTER")
  }else{
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
      colnames(info) <- c("QUAL","FILTER")
    } else {
      info_tmp_2 <- as.data.frame(
        do.call(rbind, lapply(info_list, function(rec) unname(rec[keys_all]))),
        stringsAsFactors = FALSE)
      info <- cbind(info_tmp_1, info_tmp_2)
      colnames(info) <- c("QUAL","FILTER",keys_all)
    }
  }
  
  # identify which SNPs have more than 2 alleles
  if("AC" %in% colnames(info)){
    more_alleles <- grep(",",info$AC)
    if(length(more_alleles)!=0){
      info <- info[-more_alleles,]
      # Numerify only the columns that parse as numeric; coercing every
      # column destroyed FILTER ("PASS" -> NA) and any character INFO
      # field whenever a multi-allelic record was present (F4).
      numerifiable <- vapply(info, function(col) {
        col <- as.character(col)
        all(is.na(col) | !is.na(suppressWarnings(as.numeric(col))))
      }, logical(1))
      info[numerifiable] <- lapply(info[numerifiable],
                                   function(col) as.numeric(as.character(col)))
      x@loc.all <- loc.all[-more_alleles]
      x@chromosome <- as.factor(chrom[-more_alleles])
      x@position <- pos[-more_alleles]
    }else{
      x@loc.all <- loc.all
      x@chromosome <- as.factor(chrom)
      x@position <- pos
    }
  }else{
    ALT <- vcfR::getALT(vcf)
    more_alleles <- grep(pattern = ",",ALT)
    if(length(more_alleles)>0){
      info <- info[-more_alleles,]
      x@loc.all <- loc.all[-more_alleles]
      x@chromosome <- as.factor(chrom[-more_alleles])
      x@position <- pos[-more_alleles]
    }else{
      x@loc.all <- loc.all
      x@chromosome <- as.factor(chrom)
      x@position <- pos
    }
  }

  # In genotype mode heterozygous calls are collapsed onto a 0/1/2 scale,
  # so the object is coded diploid; in dosage mode the ploidy is set from
  # the data's maximum copy number, per the documented dosage semantics --
  # forcing 2 stamped polyploid dosages (up to ploidy n) as diploid (F1).
  # The stamp is uniform across individuals because downstream checks
  # (gl.compliance.check) require a single ploidy level per object.
  if (mode == "genotype") {
    ploidy(x) <- rep(2,nInd(x))
  } else {
    ploidy(x) <- rep(max(ploidy(x)),nInd(x))
  }

  x <- gl.compliance.check(x, verbose = verbose)
  
  x$other$loc.metrics <- cbind(x$other$loc.metrics,info)
  x$other$loc.metrics$QUAL <- as.numeric(x$other$loc.metrics$QUAL)
  
  # additional metadata and long lat to the data file are stored in other
  
  if (!is.null(ind.metafile)) {
    if (verbose >= 2) {
      cat(report(
        paste("Adding individual metrics:", ind.metafile, ".\n")
      ))
    }
    ###### population and individual file to link numbers to populations...
    ind.cov <- read.csv(ind.metafile,  
                        header = TRUE, 
                        stringsAsFactors = TRUE)
    # is there an entry for every individual
    
    id.col <- match("id", names(ind.cov))
    
    if (is.na(id.col)) {
      stop(error("Fatal Error: There is no id column\n"))
    } else {
      ind.cov[, id.col] <-
        trimws(ind.cov[, id.col], which = "both")  #trim spaces
      
      if (length(ind.cov[, id.col]) != length(unique(ind.cov[, id.col]))) {
        cat(error(
          "Individual names are not unique. You need to change them!\n"
        ))
        stop()
      }
      
      # reorder
      # The id-mismatch listing affects the result, so it prints from
      # verbose >= 1 rather than unconditionally (F3).
      if (length(ind.cov[, id.col]) != length(indNames(x)) & verbose >= 1) {
        cat(
          warn(
            "Ids for individual metadata does not match the number of ids in the SNP data file. Maybe this is fine if a subset matches.\n"
          )
        )
        nam.indmeta <- ind.cov[, id.col]
        nam.dart <- indNames(x)

        nm.indmeta <- nam.indmeta[!nam.indmeta %in% nam.dart]
        nm.inddart <- nam.dart[!nam.dart %in% nam.indmeta]
        if (length(nm.indmeta) > 0) {
          cat(warn("ind.metafile ids not matched were:\n"))
          print(nm.indmeta)
        }
        if (length(nm.inddart) > 0) {
          cat(warn("DArT file ids not matched were:\n"))
          print(nm.inddart)
        }
      }

      # Metadata rows are aligned to the individuals in the SNP data by
      # id; ord carries NA for individuals without a metafile entry.
      ord <- match(indNames(x), ind.cov[, id.col])

      if (sum(!is.na(ord)) > 1) {
        if (verbose >= 2) {
          cat(report(
            paste(
              "  Ids for individual metadata (at least a subset of) are matching!\n"
            )
          ))
          cat(report(
            paste(
              "  Found ",
              sum(!is.na(ord)),
              "matching ids out of",
              nrow(ind.cov),
              "ids provided in the ind.metadata file.\n "
            )
          ))
        }
        # Individuals absent from the metafile are retained with NA
        # metadata rather than silently dropped from the object (F5).
        if (any(is.na(ord)) & verbose >= 1) {
          cat(warn(
            "  Warning: individuals without an entry in the ind.metafile are retained with NA metadata:\n"
          ))
          print(indNames(x)[is.na(ord)])
        }
      } else {
        stop(error(
          "Fatal Error: Individual ids are not matching!!!!\n"
        ))
      }
    }
    
    pop.col <- match("pop", names(ind.cov))
    
    if (is.na(pop.col)) {
      if (verbose >= 1) {
        cat(
          warn(
            "  Warning: There is no pop column, created one with all pop1 as default for all individuals\n"
          )
        )
      }
      pop(x) <- factor(rep("pop1", nInd(x)))
    } else {
      pop(x) <- as.factor(ind.cov[ord, pop.col])
      if (verbose >= 2) {
        cat(report(" Added population assignments.\n"))
      }
    }
    
    lat.col <- match("lat", names(ind.cov))
    lon.col <- match("lon", names(ind.cov))
    if (verbose >= 2) {
      if (is.na(lat.col)) {
        cat(
          warn(
            "Warning: Individual metrics do not include a latitude (lat) column\n"
          )
        )
      }
      if (is.na(lon.col)) {
        cat(
          warn(
            "Warning: Individual metrics do not include a longitude (lon) column\n"
          )
        )
      }
    }
    if (!is.na(lat.col) & !is.na(lon.col)) {
      x@other$latlon <- ind.cov[ord, c(lat.col, lon.col)]
      # Row names come from the object, not the metafile, so retained
      # individuals without a metafile entry (NA in ord) keep valid
      # row names (F5).
      rownames(x@other$latlon) <- indNames(x)
      if (verbose >= 2) {
        cat(report("  Added latlon data.\n"))
      }
    }

    other.col <- names(ind.cov)
    if (length(other.col) > 0) {
      x@other$ind.metrics <- ind.cov[ord, other.col, drop = FALSE]
      # id filled from the object so retained individuals without a
      # metafile entry carry their id rather than NA (F5).
      x@other$ind.metrics[, id.col] <- indNames(x)
      rownames(x@other$ind.metrics) <- indNames(x)
      if (verbose >= 2) {
        cat(report(
          paste(
            " Added ",
            other.col,
            " to the other$ind.metrics slot.\n"
          )
        ))
      }
    }
  }
  
  # add history
  x@other$history <- list(match.call())
  x <- gl.recalc.metrics(x, verbose = verbose)
  
  if (verbose > 2) {
    cat(
      important(
        "Genlight object does not have individual metrics. You need to add them 'manually' to the @other$ind.metrics slot.\n"
      )
    )
  }
  #convert to fbm 
  if (fbm) {
  x <- gl.gen2fbm(x, verbose = verbose) 
  }
  # else x@fbm <- NULL
  if (verbose>2) {
    cat(report(" Created an  file-backed matrix (fbm) dartR object\n"))
  } 
  
  
  return(x)
  
}
