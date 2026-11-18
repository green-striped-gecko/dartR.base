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
#' @param read.counts Logical, if TRUE the reference and alternative read
#' counts of each genotype are stored in the dartR object (see details). They
#' take much more memory than the genotypes [default FALSE].
#' @param likelihoods Logical, if TRUE the likelihoods of the three genotypes
#' (FORMAT PL, or GL) of each genotype call are stored in the dartR object
#' (see details). Only in "genotype" mode. They take 12 bytes per genotype,
#' about 190 MB for 171 individuals and 93,000 loci [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @param parser "vcfR" reads the vcf with package vcfR (the whole file is
#' held in memory as text); "cpp" reads the genotypes with the C++ parser of
#' package tidypopgen, chunk by chunk straight into a file-backed matrix,
#' which is much faster and needs little memory for large files. tidypopgen
#' is not installed with dartR.base; install it with
#' \code{install.packages("tidypopgen")} to use "cpp". See details for the
#' limits of "cpp" [default "vcfR"].
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
#'
#' When the vcf holds FORMAT DP, the mean read depth of each locus over its
#' called genotypes is stored in \code{x@other$loc.metrics$rdepth}, where
#' \code{\link{gl.report.rdepth}} and \code{\link{gl.filter.rdepth}} find it.
#' Missing and half-missing calls (e.g. ./., 0/.) are not counted; some
#' callers report a depth of 0 for them.
#'
#' With \code{read.counts = TRUE}, \code{x@other$geno.metrics} holds two
#' integer matrices of individuals x loci: \code{ref.count} and
#' \code{alt.count}, the reads of the REF and ALT allele of each genotype.
#' They are read from FORMAT AD or, for genotypes without AD, from FORMAT
#' BSDP (read counts of the bases A, C, G and T) at the REF and ALT bases.
#' Counts are stored as the vcf reports them, also for uncalled genotypes.
#'
#' With \code{likelihoods = TRUE}, \code{x@other$geno.metrics} also holds
#' \code{pl.ref}, \code{pl.het} and \code{pl.alt}: the phred-scaled
#' likelihoods of 0/0, 0/1 and 1/1 for each genotype, as integer matrices of
#' individuals x loci. They are read from FORMAT PL or, for genotypes without
#' PL, from FORMAT GL as round(-10 * GL). Each genotype's likelihoods are
#' shifted so that the most likely genotype has 0, as some callers write PL
#' whose lowest value is above 0; likelihoods are defined only up to a
#' constant, so nothing is lost. A likelihood of 0 (GL of -inf) is stored as
#' the largest integer. Genotypes whose likelihoods are missing, do not
#' number three (haploid or polyploid calls) or are all equal (e.g.
#' ./.:0,0,0) get NA; uncalled genotypes with informative likelihoods keep
#' them. Likelihoods are not stored in "dosage" mode.
#'
#' The matrices in \code{x@other$geno.metrics} have no row or column names:
#' row i belongs to \code{indNames(x)[i]} and column j to
#' \code{locNames(x)[j]}, so renaming individuals or loci cannot leave stale
#' names on them. The matrices follow the genotypes when individuals or loci
#' are subset, sorted or combined, and \code{\link{gl.randomize.snps}} swaps
#' the REF and ALT read counts and the 0/0 and 1/1 likelihoods of the loci it
#' recodes; \code{\link{gl.compliance.check}} removes them if their
#' dimensions no longer match the object.
#'
#' With \code{parser = "cpp"} (package tidypopgen) the genotypes are read
#' by tidypopgen's C++ parser into a file-backed matrix and copied in blocks
#' of loci into the returned object (file-backed with \code{fbm = TRUE},
#' converted to @gen in blocks otherwise); the fixed columns of the vcf
#' (CHROM to INFO) are read in chunks as well, so the whole file is never
#' held in memory. For the variants it reads, the result is the same as with
#' \code{parser = "vcfR"}: genotypes (count of the ALT allele), locus names,
#' alleles, chromosome, position and the QUAL, FILTER and INFO locus
#' metrics. Limits: only single-base (SNP) variants are read, so indels are
#' skipped (with a warning); only \code{mode = "genotype"}; no rdepth, read
#' counts or likelihoods; and only diploid genotype calls -- a vcf with
#' haploid or polyploid calls (e.g. a sex chromosome in the heterogametic
#' sex) stops with an error, as tidypopgen would code a haploid call like a
#' heterozygote. Use \code{parser = "vcfR"} in these cases.
#' @return A genlight object.
#' @export
#' @author Author(s): Bernd Gruber, Ching Ching Lau. Custodian: Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
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
                        read.counts = FALSE,
                        likelihoods = FALSE,
                        verbose = NULL,
                        parser = "vcfR") {

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  funname <- utils.flag.start(func = funname,
                   build = "Jackson",
                   verbose = verbose)
  
  x <- NULL

  if (length(parser) != 1 || !parser %in% c("vcfR", "cpp")) {
    stop(error("  Fatal Error: parser must be \"vcfR\" or \"cpp\"\n"))
  }
  if (parser == "cpp") {
    if (!requireNamespace("tidypopgen", quietly = TRUE)) {
      stop(error("  Fatal Error: parser = \"cpp\" needs the package",
                 "tidypopgen, which is not installed with dartR.base.",
                 "Install it with install.packages(\"tidypopgen\"),",
                 "or use parser = \"vcfR\"\n"))
    }
    if (mode != "genotype" || read.counts || likelihoods) {
      stop(error("  Fatal Error: parser = \"cpp\" reads genotypes only",
                 "(mode = \"genotype\", no read counts or likelihoods);",
                 "use parser = \"vcfR\" for these options\n"))
    }
    # genotypes by tidypopgen into a file-backed matrix, fixed columns in
    # chunks (see .vcf_read_cpp); no read depth from this route
    rc <- .vcf_read_cpp(vcffile, verbose = verbose)
    x <- rc$x
    info <- rc$info
    rm(rc)
    depth <- list()
  } else {
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
  
  # PL holds one likelihood per copy number of the ALT allele, so only
  # the three diploid genotypes of genotype mode are kept
  if (likelihoods && mode != "genotype") {
    if (verbose >= 1) {
      cat(warn("  Warning: genotype likelihoods are stored only in 'genotype' mode, so none were stored.\n"))
    }
    likelihoods <- FALSE
  }

  # Depth, read counts and likelihoods are read from the same biallelic
  # records that became the loci, so they line up with the genotypes
  depth <- utils.vcf.depth(vcf, loci = which(vcfR::is.biallelic(vcf)),
                           read.counts = read.counts,
                           likelihoods = likelihoods)

  # adding SNP information from VCF (QUAL, FILTER and INFO of the
  # biallelic records; shared with the cpp parser)
  fi <- .vcf_fix_info(vcf@fix, pos = pos, loc.all = loc.all)
  info <- fi$info
  x@loc.all <- fi$loc.all
  x@chromosome <- fi$chromosome
  x@position <- fi$position
  rm(vcf, fi)
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
  
  if (!is.null(depth$rdepth)) {
    x$other$loc.metrics$rdepth <- depth$rdepth
    if (verbose >= 2) {
      cat(report("  Added rdepth (mean FORMAT DP over called genotypes) to the other$loc.metrics slot.\n"))
    }
  }
  if (read.counts) {
    if (is.null(depth$ref.count)) {
      if (verbose >= 1) {
        cat(warn("  Warning: the vcf has no AD or BSDP field, so no read counts were stored.\n"))
      }
    } else {
      x@other$geno.metrics <- list(ref.count = depth$ref.count,
                                   alt.count = depth$alt.count)
      if (verbose >= 2) {
        cat(report("  Added REF and ALT read counts (ref.count, alt.count) to the other$geno.metrics slot.\n"))
      }
    }
  }
  if (likelihoods) {
    if (is.null(depth$pl.ref)) {
      if (verbose >= 1) {
        cat(warn("  Warning: the vcf has no PL or GL field, so no genotype likelihoods were stored.\n"))
      }
    } else {
      x@other$geno.metrics <- c(x@other$geno.metrics,
                                depth[c("pl.ref", "pl.het", "pl.alt")])
      if (verbose >= 2) {
        cat(report("  Added genotype likelihoods (pl.ref, pl.het, pl.alt) to the other$geno.metrics slot.\n"))
      }
    }
  }
  rm(depth)

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
  x@other$history <- list(utils.history.call(match.call()))
  x <- gl.recalc.metrics(x, verbose = verbose)
  
  if (verbose > 2) {
    cat(
      important(
        "Genlight object does not have individual metrics. You need to add them 'manually' to the @other$ind.metrics slot.\n"
      )
    )
  }
  # convert to fbm: the cpp parser reads into a file-backed matrix, which
  # is kept with fbm = TRUE and converted in blocks of loci otherwise
  if (fbm) {
    if (!.has_fbm(x)) x <- gl.gen2fbm(x, verbose = verbose)
  } else if (.has_fbm(x)) {
    x <- gl.fbm2gen(x, verbose = 0)
  }
  # else x@fbm <- NULL
  if (verbose>2) {
    cat(report(" Created an  file-backed matrix (fbm) dartR object\n"))
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(x)
  
}
