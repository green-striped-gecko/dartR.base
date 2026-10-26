#' @name utils.vcf.depth
#' @title Utility function to read locus depth, allele read counts and
#' genotype likelihoods from a vcfR object
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY
#' END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.
#' @param vcf A vcfR object [defined in function \code{\link{gl.read.vcf}}].
#' @param loci Indices of the records kept as loci, in the order of the loci
#' in the genlight object [default all records].
#' @param read.counts If TRUE, the reference and alternative read counts of
#' each genotype are also returned [default FALSE].
#' @param likelihoods If TRUE, the likelihoods of the three diploid genotypes
#' of each genotype call are also returned [default FALSE].
#' @param chunk Number of records processed at a time, which bounds the
#' memory used [default 10000].
#' @details
#' rdepth is the mean of FORMAT DP over the called genotypes of a locus.
#' Missing and half-missing calls (e.g. ./., 0/.) are not called, as in
#' \code{\link{utils.vcfr2genlight.polyploid}}; some callers report DP = 0
#' for them. A locus without called genotypes or without DP gets NA.
#'
#' Read counts come from FORMAT AD (first value REF, second ALT). Where a
#' genotype has no AD, they come from FORMAT BSDP (read counts of the bases A,
#' C, G and T, as written by some callers) taken at the REF and ALT bases;
#' records whose REF or ALT is not a single base get NA. Counts are stored as
#' reported, also for uncalled genotypes.
#'
#' Likelihoods come from FORMAT PL or, for genotypes without PL, from FORMAT
#' GL as round(-10 * GL). They are phred-scaled and shifted so that the most
#' likely genotype has 0; a likelihood of 0 (-inf in GL) becomes the largest
#' integer. A genotype whose likelihoods are missing, do not number three or
#' are all equal (e.g. ./.:0,0,0) gets NA; uncalled genotypes with
#' informative likelihoods keep them.
#' @return A list with rdepth (one value per locus, NULL when the vcf has no
#' DP); when read.counts is TRUE, ref.count and alt.count (integer matrices,
#' individuals x loci, NULL when the vcf has neither AD nor BSDP); and when
#' likelihoods is TRUE, pl.ref, pl.het and pl.alt (the phred-scaled
#' likelihoods of 0/0, 0/1 and 1/1, integer matrices of individuals x loci,
#' NULL when the vcf has neither PL nor GL).
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' depth <- utils.vcf.depth(vcf, read.counts = TRUE, likelihoods = TRUE)
#' }
#' @export

utils.vcf.depth <- function(vcf,
                            loci = seq_len(nrow(vcf@fix)),
                            read.counts = FALSE,
                            likelihoods = FALSE,
                            chunk = 10000L) {
  res <- list(rdepth = NULL, ref.count = NULL, alt.count = NULL,
              pl.ref = NULL, pl.het = NULL, pl.alt = NULL)
  if (ncol(vcf@gt) < 2 || length(loci) == 0) {
    return(res)
  }

  fields <- unique(unlist(strsplit(unique(vcf@gt[loci, "FORMAT"]), ":",
                                   fixed = TRUE)))
  has.dp <- "DP" %in% fields
  has.ad <- read.counts && "AD" %in% fields
  has.bsdp <- read.counts && "BSDP" %in% fields
  has.pl <- likelihoods && "PL" %in% fields
  has.gl <- likelihoods && "GL" %in% fields
  n.ind <- ncol(vcf@gt) - 1L
  n.loc <- length(loci)

  if (has.dp) {
    res$rdepth <- rep(NA_real_, n.loc)
  }
  if (has.ad || has.bsdp) {
    res$ref.count <- matrix(NA_integer_, n.ind, n.loc)
    res$alt.count <- matrix(NA_integer_, n.ind, n.loc)
  }
  pl.names <- c("pl.ref", "pl.het", "pl.alt")
  if (has.pl || has.gl) {
    res[pl.names] <- list(matrix(NA_integer_, n.ind, n.loc))
  }
  if (!has.dp && is.null(res$ref.count) && is.null(res$pl.ref)) {
    return(res)
  }

  bases <- c("A", "C", "G", "T")
  for (s in seq(1L, n.loc, by = chunk)) {
    cols <- s:min(n.loc, s + chunk - 1L)
    v <- vcf[loci[cols], ]

    if (has.dp) {
      gt <- vcfR::extract.gt(v, element = "GT", IDtoRowNames = FALSE)
      dp <- vcfR::extract.gt(v, element = "DP", as.numeric = TRUE,
                             IDtoRowNames = FALSE)
      dp[is.na(gt) | grepl(".", gt, fixed = TRUE)] <- NA
      n.called <- rowSums(!is.na(dp))
      res$rdepth[cols] <- ifelse(n.called > 0,
                                 rowSums(dp, na.rm = TRUE) / n.called, NA)
    }

    if (!is.null(res$ref.count)) {
      ref <- alt <- matrix(NA_real_, length(cols), n.ind)
      if (has.ad) {
        ad <- vcfR::extract.gt(v, element = "AD", IDtoRowNames = FALSE)
        ref <- vcfR::masplit(ad, record = 1, sort = 0)
        alt <- vcfR::masplit(ad, record = 2, sort = 0)
      }
      # genotypes without AD take the counts of their REF and ALT bases
      need <- is.na(ref) & is.na(alt)
      if (has.bsdp && any(need)) {
        bsdp <- vcfR::extract.gt(v, element = "BSDP", IDtoRowNames = FALSE)
        ref.base <- match(v@fix[, "REF"], bases)
        alt.base <- match(v@fix[, "ALT"], bases)
        # at an indel the base counts do not separate the two alleles, even
        # when one allele is a single base
        snv <- !is.na(ref.base) & !is.na(alt.base)
        ref.base[!snv] <- NA
        alt.base[!snv] <- NA
        for (k in seq_along(bases)) {
          if (!any(ref.base %in% k | alt.base %in% k)) next
          count <- vcfR::masplit(bsdp, record = k, sort = 0)
          # a vector of one value per record recycles down the columns
          sel <- need & (ref.base %in% k)
          ref[sel] <- count[sel]
          sel <- need & (alt.base %in% k)
          alt[sel] <- count[sel]
        }
      }
      res$ref.count[, cols] <- as.integer(t(ref))
      res$alt.count[, cols] <- as.integer(t(alt))
    }

    if (!is.null(res$pl.ref)) {
      lik <- rep(list(matrix(NA_real_, length(cols), n.ind)), 3L)
      if (has.pl) {
        pl <- vcfR::extract.gt(v, element = "PL", IDtoRowNames = FALSE)
        # "." and single values cannot hold three likelihoods, and masplit
        # reads an empty field as 0
        pl[!grepl(",", pl, fixed = TRUE)] <- NA
        three <- vcfR::masplit(pl, count = TRUE, sort = 0) %in% 3
        for (k in 1:3) {
          lik[[k]] <- vcfR::masplit(pl, record = k, sort = 0)
          lik[[k]][!three] <- NA
        }
      }
      # genotypes without PL take their GL; it is split here because
      # masplit reads -inf as 0
      need <- is.na(lik[[1]]) | is.na(lik[[2]]) | is.na(lik[[3]])
      if (has.gl && any(need)) {
        gl <- vcfR::extract.gt(v, element = "GL", IDtoRowNames = FALSE)
        sel <- need & !is.na(gl)
        parts <- strsplit(gl[sel], ",", fixed = TRUE)
        three <- lengths(parts) == 3L
        vals <- matrix(NA_real_, length(parts), 3L)
        vals[three, ] <- matrix(suppressWarnings(
          as.numeric(unlist(parts[three]))), ncol = 3L, byrow = TRUE)
        for (k in 1:3) {
          lik[[k]][sel] <- round(-10 * vals[, k])
        }
      }
      # likelihoods are relative, so the most likely genotype is set to 0;
      # equal likelihoods (e.g. ./.:0,0,0) carry no information
      low <- do.call(pmin, lik)
      low[which(low == do.call(pmax, lik))] <- NA
      for (k in 1:3) {
        p <- pmin(lik[[k]] - low, .Machine$integer.max)
        res[[pl.names[k]]][, cols] <- as.integer(t(p))
      }
    }
  }

  res
}
