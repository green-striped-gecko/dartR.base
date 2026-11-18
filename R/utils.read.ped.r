#' @name utils.read.ped
#' @title Reads a PLINK .ped file into a SnpMatrix with family and map data
#' @family utilities
#'
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.
#'
#' @details
#' A vendored, lightly modified copy of \code{snpStats::read.pedfile}. The
#' file is read twice: a first pass counts the lines, a second parses each
#' line into the six pedigree columns and the genotype allele pairs. For
#' each locus the first allele encountered in file order becomes
#' \code{allele.1} and the second becomes \code{allele.2}; genotypes are
#' coded on that basis (raw SnpMatrix coding; as numeric, the count of
#' \code{allele.2}). Loci at which more than two alleles are observed are
#' set entirely to NA. The only in-package caller is
#' \code{gl.report.ld.map}.
#'
#' @param file Name of the .ped file to read (plain text or gzipped)
#' [required].
#' @param snps Either the name of the associated .map file, or a character
#' vector of locus names (one per locus) [optional; if missing, loci are
#' named locus.1, locus.2, ...].
#' @param which When snps is a map file, the column of that file holding
#' the locus names [optional; defaults to the first column with no
#' duplicates].
#' @param split Regular expression used to split fields on each line
#' [default "\\t| +", i.e. tabs and spaces].
#' @param sep Separator used when constructing locus and subject names
#' [default "."].
#' @param na.strings Strings to be treated as missing values. Note that
#' the default treats the conventional PLINK missing-allele code as
#' missing [default "0"].
#' @param lex.order If TRUE, alleles at each locus are reordered
#' lexicographically, and the genotype codes are switched to match
#' [default FALSE].
#' @param show_warnings If FALSE, the warnings for no-data, monomorphic
#' and multi-allelic loci are suppressed. The multi-allelic NA reset is
#' applied regardless [default TRUE].
#'
#' @return A list with elements: \code{genotypes} (a
#' \code{snpStats::SnpMatrix}, individuals x loci), \code{fam} (a data
#' frame with pedigree, member, father, mother, sex, affected) and
#' \code{map} (a data frame with the locus names and alleles, plus the
#' .map columns when one was supplied).
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @keywords internal
#' @export
#' @importFrom snpStats switch.alleles

utils.read.ped <- function (file, 
                            # n, 
                            snps, 
                            which, 
                            split = "\t| +",
                            sep = ".", 
                            na.strings = "0", 
                            lex.order = FALSE,
                            show_warnings = TRUE){
  
# check if packages are installed
  pkg <- "snpStats"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }



  r0 <- as.raw(0)
  r1 <- as.raw(1)
  r2 <- as.raw(2)
  r3 <- as.raw(3)
  con <- gzfile(file)
  open(con)
  # if (missing(n)) {
  n <- 0
  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) 
      break
    n <- n + 1
  }
  if (n == 0){ 
    stop("Nothing read")
  }
  seek(con, 0)
  # }
  gen <- missing(snps)
  map <- NULL
  if (!gen) {
    m <- length(snps)
    if (m == 1) {
      map <- read.table(snps, comment.char = "")
      m <- nrow(map)
      if (missing(which)) {
        which <- 1
        repeat {
          snps <- map[, which]
          if (!any(duplicated(snps))) 
            break
          if (which == ncol(map)) 
            stop("No unambiguous snp names found on file")
          which <- which + 1
        }
      } else {
        snps <- map[, which]
      }
    }
  } else {
    line <- readLines(con, n = 1)
    fields <- strsplit(line, split)[[1]]
    nf <- length(fields)
    if (nf%%2 != 0) 
      stop("Odd number of fields")
    m <- (nf - 6)/2
    seek(con, 0)
  }
  nf <- 6 + 2 * m
  result <- matrix(raw(n * m), nrow = n)
  ped <- character(n)
  mem <- character(n)
  pa <- character(n)
  ma <- character(n)
  sex <- numeric(n)
  aff <- numeric(n)
  rownms <- character(n)
  a1 <- a2 <- rep(NA, m)
  a1m <- a2m <- rep(TRUE, m)
  mallelic <- rep(FALSE, m)
  for (i in 1:n) {
    line <- readLines(con, n = 1)
    fields <- strsplit(line, split)[[1]]
    to.na <- fields %in% na.strings
    fields[to.na] <- NA
    ped[i] <- fields[1]
    mem[i] <- fields[2]
    pa[i] <- fields[3]
    ma[i] <- fields[4]
    sex[i] <- suppressWarnings(as.numeric(fields[5]))
    aff[i] <- as.numeric(fields[6])
    alleles <- matrix(fields[7:nf], byrow = TRUE, ncol = 2)
    one <- two <- rep(FALSE, m)
    for (k in 1:2) {
      ak <- alleles[, k]
      akm <- is.na(ak)
      br1 <- !akm & a1m
      a1[br1] <- ak[br1]
      a1m[br1] <- FALSE
      br2 <- !akm & (a1 == ak)
      one[br2] <- TRUE
      br3 <- !akm & !a1m & (a1 != ak)
      br4 <- br3 & a2m
      a2[br4] <- ak[br4]
      a2m[br4] <- FALSE
      br5 <- br3 & (a2 == ak)
      two[br5] <- TRUE
      # A locus is multi-allelic when THIS column's allele matches
      # neither a1 (br2) nor a2 (br5) and is not missing. The original
      # test used the accumulated one/two, which are set per individual
      # across both allele columns, so a novel allele paired with a
      # known one escaped detection (F3).
      mallelic <- mallelic | !(akm | br2 | br5)
    }
    gt <- rep(r0, m)
    gt[one & !two] <- r1
    gt[one & two] <- r2
    gt[two & !one] <- r3
    result[i, ] <- gt
  }
  close(con)
  if (any(a1m & show_warnings==T)){ 
    warning("no data for ", sum(a1m), " loci")
  }
  mono <- (a2m & !a1m)
  if (any(mono & show_warnings==T)){ 
    warning(sum(mono), " loci were monomorphic")
  }
  if (any(mallelic)) {
    # The NA reset is data cleaning, not messaging: it must run
    # regardless of show_warnings; only the warning is gated (F2).
    result[, mallelic] <- r0
    if (show_warnings) {
      warning(sum(mallelic), " loci were multi-allelic --- set to NA")
    }
  }
  if (gen){ 
    snps <- paste("locus", 1:m, sep = sep)
  }
  if (any(duplicated(ped))) {
    if (any(duplicated(mem))) {
      rnames <- paste(ped, mem, sep = sep)
      if (any(duplicated(rnames))){ 
        stop("could not create unique subject identifiers")
      }
    } else{
      rnames <- mem
    }
  }  else {
    rnames <- ped
  }
  dimnames(result) <- list(rnames, snps)
  result <- new("SnpMatrix", result)
  if (lex.order) {
    swa <- (!(is.na(a1) | is.na(a2)) & (a1 > a2))
    # Assign the switched matrix; the original discarded the return
    # value, swapping the map alleles but not the genotypes (F1).
    result <- snpStats::switch.alleles(result, swa)
    a1n <- a1
    a1n[swa] <- a2[swa]
    a2[swa] <- a1[swa]
    a1 <- a1n
  }
  fam <- data.frame(row.names = rnames, pedigree = ped, member = mem, 
                    father = pa, mother = ma, sex = sex, affected = aff)
  if (is.null(map)){ 
    map <- data.frame(row.names = snps, snp.name = snps, 
                      allele.1 = a1, allele.2 = a2)
  } else {
    map$allele.1 <- a1
    map$allele.2 <- a2
    names(map)[which] <- "snp.names"
  }
  list(genotypes = result, fam = fam, map = map)
}

rotate.matrix <- function (x, 
                           angle = 10, 
                           method = "bilinear"){
  
  img <- x
  angle.rad <- angle * pi/180
  co.x <- matrix(rep(-(ncol(img)/2 - 0.5):(ncol(img)/2 - 
                                             0.5), nrow(img)), nrow = nrow(img), byrow = T)
  co.y <- matrix(rep(-(nrow(img)/2 - 0.5):(nrow(img)/2 - 
                                             0.5), ncol(img)), ncol = ncol(img))
  co.xn <- round(co.x * cos(angle.rad) - co.y * sin(angle.rad))
  co.yn <- round(co.x * sin(angle.rad) + co.y * cos(angle.rad))
  co.xn2 <- co.xn + max(co.xn) + 1
  co.yn2 <- co.yn + max(co.yn) + 1
  img.rot <- numeric(max(co.yn2) * max(co.xn2))
  img.rot[(co.xn2 - 1) * max(co.yn2) + co.yn2] <- img
  dim(img.rot) <- c(max(co.yn2), max(co.xn2))
  attr(img.rot, "bits.per.sample") <- attr(img, "bits.per.sample")
  attr(img.rot, "samples.per.pixel") <- attr(img, "samples.per.pixel")
  return(img.rot)
  
}