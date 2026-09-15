#' @name gl.report.contamination
#' @title Screens a genlight object for cross-contaminated samples
#' @family unmatched report
#'
#' @description
#' Flags individuals whose genotypes carry DNA from another sample, using two
#' tiers of evidence that need only the genotypes, the population labels and,
#' when available, the plate wells stored by \code{gl.read.dart()}.
#'
#' @details
#' Tier 1 scores each individual against the other members of its own
#' population. Heterozygosity is the proportion of called loci scored
#' heterozygous. The rare-allele burden is the proportion of called loci at
#' which the individual carries an allele whose frequency among the other
#' members of its population is below \code{rare.freq}; loci with fewer than
#' \code{min.n} other individuals called are skipped. Both statistics are
#' reported as a robust z-score (median and MAD within the population) and as
#' a fold over the population median.
#'
#' Tier 2 names the most likely source of the foreign DNA. For an individual
#' with at least \code{min.share} rare-allele loci, every other individual is
#' scored by the proportion of those rare alleles it carries (sharing). The
#' candidates within \code{share.tol} of the best score, which in practice are
#' the members of one population, are then separated by residual kinship: a
#' VanRaden genomic relationship matrix halved to kinship, with each pair
#' expressed as a residual over the median kinship of its population pair so
#' that within- and between-population pairs are on the same scale. For an
#' individual with too few rare alleles (contamination from its own
#' population, or a clean sample) the partner is simply the largest residual.
#' The partner's sharing score and its residual kinship z-score against that
#' individual's other residuals are reported. If plate positions are known,
#' the partner is tested for plate adjacency (same plate, wells sharing an
#' edge). Positions are taken from
#' \code{plate} if supplied, otherwise from the \code{plate_location} column
#' of \code{@@other$ind.metrics} written by \code{gl.read.dart()} (for
#' example "1-C2"), otherwise from \code{plate} and \code{well} columns.
#'
#' Flags. An individual is a "suspect" when its heterozygosity has z-score
#' above \code{z.flag} and exceeds the population median by more than
#' \code{min.excess}.
#' Contamination always raises heterozygosity, so this is the required
#' signal. A suspect whose top kinship partner sits in an adjacent well is
#' upgraded to "adjacent". An individual whose rare-allele burden is elevated
#' by the same criteria but whose heterozygosity is normal is reported as
#' "rare-only" (z-score above \code{z.flag} and excess above
#' \code{rare.min.excess}): this is the signature of population structure
#' inside the assigned population, or of a mislabelled individual, not of
#' contamination.
#'
#' Limits. The screen is a candidate list, not a verdict. Contamination
#' from the same population as the host raises only heterozygosity and
#' kinship, so a low rare-allele burden does not clear a suspect, and a
#' population whose members differ widely in heterozygosity (inbreeding,
#' admixture, relatives) can yield suspects that are not contaminated.
#' Power falls as host and contaminant become genetically closer: a mixture
#' between two populations of one species at a typical heterozygosity of
#' 0.3 shifts heterozygosity by less than the natural spread of many
#' populations. The decisive test in that case is allele balance from read
#' counts, which this function does not use. Tier 2 can only name a source
#' that is in the dataset.
#'
#' The function densifies the genotype matrix with \code{as.matrix()}.
#'
#' The plot shows heterozygosity against rare-allele burden with suspects
#' labelled.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param rare.freq Allele frequency in the rest of the population below
#' which an allele counts as rare [default 0.02].
#' @param min.n Minimum number of other individuals in the population called
#' at a locus for that locus to enter the rare-allele burden [default 5].
#' @param z.flag Robust z-score above which a statistic is an outlier
#' [default 3].
#' @param min.excess Absolute amount by which heterozygosity must also exceed
#' the population median to be flagged; a robust z-score alone over-flags
#' populations whose members are nearly identical [default 0.02].
#' @param rare.min.excess The same absolute margin for the rare-allele
#' burden, which is an order of magnitude smaller than heterozygosity
#' [default 0.005].
#' @param min.share Minimum number of rare-allele loci an individual needs
#' before its source is chosen by allele sharing rather than by residual
#' kinship alone [default 20].
#' @param share.tol Sharing scores within this distance of the best are
#' treated as ties and separated by residual kinship [default 0.05].
#' @param plate Data frame with columns id, plate and well (for example
#' "C4") giving plate positions; overrides positions found in the individual
#' metadata [default NULL].
#' @param plot.display If TRUE, resultant plots are displayed in the plot
#' window [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors List of two color names, the first for suspect
#' individuals and the second for the rest [default c("#2171B5","#6BAED6")].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#' @param ... Parameters passed to function \link[ggplot2]{ggsave},
#' such as width and height, when the ggplot is to be saved.
#'
#' @return A list with three elements: \code{ind}, a data frame with one row
#' per individual holding the tier 1 and tier 2 statistics (including
#' \code{n.foreign}, the number of rare-allele loci, and \code{share}, the
#' partner's sharing score) and the flag;
#' \code{pairs}, a data frame of plate-adjacent pairs with their kinship and
#' residual kinship, or NULL when no plate positions are known; and
#' \code{kinship}, the matrix of residual kinship.
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' out <- gl.report.contamination(testset.gl)
#' head(out$ind)
#'
#' @seealso \code{\link{gl.report.heterozygosity}},
#' \code{\link{gl.filter.heterozygosity}}
#' @importFrom stats mad median sd ave
#' @importFrom methods is
#' @export

gl.report.contamination <- function(x,
                                    rare.freq = 0.02,
                                    min.n = 5,
                                    z.flag = 3,
                                    min.excess = 0.02,
                                    rare.min.excess = 0.005,
                                    min.share = 20,
                                    share.tol = 0.05,
                                    plate = NULL,
                                    plot.display = TRUE,
                                    plot.theme = theme_dartR(),
                                    plot.colors = NULL,
                                    plot.dir = NULL,
                                    plot.file = NULL,
                                    verbose = NULL,
                                    ...) {

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  if (verbose == 0) plot.display <- FALSE

  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)

  # SET COLOURS
  if (is.null(plot.colors)) {
    plot.colors <- c("#2171B5", "#6BAED6")
  } else if (length(plot.colors) > 2) {
    if (verbose >= 2) {
      cat(warn("  More than 2 colors specified, only the first 2 are used\n"))
    }
    plot.colors <- plot.colors[1:2]
  }

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  if (datatype != "SNP") {
    stop(error("  Heterozygosity is undefined for presence/absence data; ",
               "this function needs SNP data\n"))
  }
  if (!is(x, "dartR")) {
    class(x) <- "dartR"
    if (verbose > 2) {
      cat(warn("  Warning: Standard adegenet genlight object encountered. ",
               "Converted to compatible dartR genlight object\n"))
    }
  }

  # FUNCTION SPECIFIC ERROR CHECKING
  if (is.null(pop(x)) || nPop(x) == 0) {
    stop(error("  Populations must be assigned before screening; ",
               "see gl.define.pop() or gl.reassign.pop()\n"))
  }
  if (rare.freq <= 0 || rare.freq >= 0.5) {
    if (verbose >= 2) {
      cat(warn("  Warning: rare.freq must lie in (0, 0.5), set to 0.02\n"))
    }
    rare.freq <- 0.02
  }
  if (min.n < 2) {
    if (verbose >= 2) {
      cat(warn("  Warning: min.n must be at least 2, set to 2\n"))
    }
    min.n <- 2
  }
  if (z.flag <= 0 || min.excess < 0 || rare.min.excess < 0) {
    stop(error("  z.flag must be positive and the excess margins ",
               "non-negative\n"))
  }
  if (min.share < 1 || share.tol < 0 || share.tol > 1) {
    stop(error("  min.share must be at least 1 and share.tol in [0, 1]\n"))
  }
  if (!is.null(plate) && !all(c("id", "plate", "well") %in% names(plate))) {
    stop(error("  plate must have columns id, plate and well\n"))
  }
  small <- names(which(table(pop(x)) < 2))
  if (length(small) > 0 && verbose >= 2) {
    cat(warn("  Warning: populations with a single individual are not ",
             "screened:", paste(small, collapse = ", "), "\n"))
  }

  # DO THE JOB
  gm <- as.matrix(x)
  ids <- indNames(x)
  pops <- as.character(pop(x))
  n.ind <- nrow(gm)
  called <- !is.na(gm)

  # Robust z-score; falls back to the SD when the MAD is zero
  robust.z <- function(v) {
    m <- median(v, na.rm = TRUE)
    s <- mad(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) s <- sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(v)))
    (v - m) / s
  }
  pop.median <- function(v) {
    ave(v, pops, FUN = function(w) median(w, na.rm = TRUE))
  }

  # Tier 1: heterozygosity and leave-one-out rare-allele burden
  if (verbose >= 2) {
    cat(report("  Tier 1: heterozygosity and rare-allele burden\n"))
  }
  het <- rowMeans(gm == 1, na.rm = TRUE)
  callrate <- rowMeans(called)
  rare.burden <- rep(NA_real_, n.ind)
  # Per individual: its rare-allele loci and the sign of the rare allele
  rare.loci <- vector("list", n.ind)
  for (p in unique(pops)) {
    idx <- which(pops == p)
    if (length(idx) < 2) next
    gp <- gm[idx, , drop = FALSE]
    alt.sum <- colSums(gp, na.rm = TRUE)
    n.call <- colSums(!is.na(gp))
    for (k in seq_along(idx)) {
      g <- gp[k, ]
      ok <- !is.na(g)
      # Allele counts among the other members of the population
      n.other <- n.call - ok
      alt.other <- alt.sum - ifelse(ok, g, 0)
      use <- ok & n.other >= min.n
      if (!any(use)) next
      p.alt <- alt.other[use] / (2 * n.other[use])
      gg <- g[use]
      rare.alt <- p.alt < rare.freq & gg >= 1
      rare.ref <- p.alt > 1 - rare.freq & gg <= 1
      rare.burden[idx[k]] <- mean(rare.alt | rare.ref)
      loc.use <- which(use)
      rare.loci[[idx[k]]] <- list(
        loci = c(loc.use[rare.alt], loc.use[rare.ref]),
        sign = c(rep(1, sum(rare.alt)), rep(-1, sum(rare.ref))))
    }
  }
  het.z <- ave(het, pops, FUN = robust.z)
  rare.z <- ave(rare.burden, pops, FUN = robust.z)
  het.excess <- het - pop.median(het)
  rare.excess <- rare.burden - pop.median(rare.burden)

  # Tier 2: residual kinship over the population-pair median
  if (verbose >= 2) {
    cat(report("  Tier 2: residual kinship and plate adjacency\n"))
  }
  col.mean <- colMeans(gm, na.rm = TRUE)
  keep <- is.finite(col.mean) & col.mean > 0 & col.mean < 2
  gi <- gm[, keep, drop = FALSE]
  gi[!called[, keep]] <- col.mean[keep][col(gi)][!called[, keep]]
  p.k <- col.mean[keep] / 2
  z.mat <- sweep(gi, 2, 2 * p.k)
  kin <- (z.mat %*% t(z.mat)) / (4 * sum(p.k * (1 - p.k)))
  dimnames(kin) <- list(ids, ids)
  resid <- kin
  for (a in unique(pops)) {
    for (b in unique(pops)) {
      ia <- pops == a
      ib <- pops == b
      blk <- kin[ia, ib, drop = FALSE]
      if (a == b) diag(blk) <- NA
      resid[ia, ib] <- kin[ia, ib] - median(blk, na.rm = TRUE)
    }
  }
  diag(resid) <- NA

  # Source: population by sharing of the individual's rare alleles, then the
  # individual within it by residual kinship
  top.j <- apply(resid, 1, which.max)
  share.top <- rep(NA_real_, n.ind)
  n.foreign <- vapply(rare.loci, function(r) length(r$loci), integer(1))
  for (i in which(n.foreign >= min.share)) {
    loci <- rare.loci[[i]]$loci
    sg <- rare.loci[[i]]$sign
    # dosage of the rare allele in every individual: gm when the rare allele
    # is the alternative, 2 - gm when it is the reference
    dos <- sweep(gm[, loci, drop = FALSE] - 1, 2, sg, `*`) + 1
    score <- rowMeans(dos, na.rm = TRUE) / 2
    score[i] <- NA
    if (all(is.na(score))) next
    cand <- which(score >= max(score, na.rm = TRUE) - share.tol)
    top.j[i] <- cand[which.max(resid[i, cand])]
    share.top[i] <- score[top.j[i]]
  }
  top.resid <- resid[cbind(seq_len(n.ind), top.j)]
  top.z <- vapply(seq_len(n.ind), function(i) {
    r <- resid[i, ]
    (top.resid[i] - median(r, na.rm = TRUE)) / mad(r, na.rm = TRUE)
  }, numeric(1))

  # Plate positions and adjacency
  pos <- NULL
  if (!is.null(plate)) {
    pos <- plate[match(ids, plate$id), c("plate", "well")]
  } else if (!is.null(x@other$ind.metrics)) {
    im <- x@other$ind.metrics
    if ("plate_location" %in% names(im)) {
      pw <- strsplit(as.character(im$plate_location), "-")
      pos <- data.frame(plate = vapply(pw, `[`, "", 1),
                        well = vapply(pw, `[`, "", 2))
    } else if (all(c("plate", "well") %in% names(im))) {
      pos <- data.frame(plate = im$plate, well = im$well)
    }
  }
  adjacent <- rep(NA, n.ind)
  pairs <- NULL
  if (!is.null(pos) && any(!is.na(pos$well))) {
    w.row <- match(substr(pos$well, 1, 1), LETTERS)
    w.col <- suppressWarnings(as.integer(substring(pos$well, 2)))
    is.adj <- function(i, j) {
      !is.na(w.row[i]) && !is.na(w.row[j]) &&
        identical(pos$plate[i], pos$plate[j]) &&
        (abs(w.row[i] - w.row[j]) + abs(w.col[i] - w.col[j])) == 1
    }
    adjacent <- vapply(seq_len(n.ind), function(i) is.adj(i, top.j[i]),
                       logical(1))
    pr <- which(upper.tri(resid), arr.ind = TRUE)
    pr <- pr[apply(pr, 1, function(q) is.adj(q[1], q[2])), , drop = FALSE]
    if (nrow(pr) > 0) {
      pairs <- data.frame(id1 = ids[pr[, 1]], id2 = ids[pr[, 2]],
                          well1 = pos$well[pr[, 1]], well2 = pos$well[pr[, 2]],
                          pop1 = pops[pr[, 1]], pop2 = pops[pr[, 2]],
                          kinship = round(kin[pr], 4),
                          resid = round(resid[pr], 4),
                          stringsAsFactors = FALSE)
      pairs <- pairs[order(-pairs$resid), ]
      rownames(pairs) <- NULL
    }
  } else if (verbose >= 2) {
    cat(warn("  No plate positions found; adjacency not tested\n"))
  }

  # Flags: heterozygosity is the required signal
  het.hit <- (het.z > z.flag & het.excess > min.excess) %in% TRUE
  rare.hit <- (rare.z > z.flag & rare.excess > rare.min.excess) %in% TRUE
  flag <- ifelse(het.hit, "suspect", ifelse(rare.hit, "rare-only", ""))
  flag[flag == "suspect" & adjacent %in% TRUE] <- "adjacent"

  ind <- data.frame(id = ids,
                    pop = pops,
                    well = if (is.null(pos)) NA else pos$well,
                    callrate = round(callrate, 3),
                    het = round(het, 4),
                    het.z = round(het.z, 1),
                    het.excess = round(het.excess, 4),
                    rare.burden = round(rare.burden, 4),
                    rare.z = round(rare.z, 1),
                    rare.excess = round(rare.excess, 4),
                    n.foreign = n.foreign,
                    top.partner = ids[top.j],
                    partner.pop = pops[top.j],
                    share = round(share.top, 3),
                    partner.well = if (is.null(pos)) NA else pos$well[top.j],
                    kin.resid = round(top.resid, 4),
                    kin.z = round(top.z, 1),
                    adjacent = adjacent,
                    flag = flag,
                    stringsAsFactors = FALSE)
  ind <- ind[order(match(ind$flag, c("adjacent", "suspect", "rare-only", "")),
                   -ind$het.excess), ]
  rownames(ind) <- NULL

  # Report
  n.sus <- sum(ind$flag %in% c("suspect", "adjacent"))
  n.rare <- sum(ind$flag == "rare-only")
  if (verbose >= 1) {
    cat(important("  Suspect individuals:", n.sus, "of", n.ind,
                  "(", sum(ind$flag == "adjacent"),
                  "with an adjacent-well partner );",
                  "rare-only:", n.rare, "\n"))
    if (n.sus > 0) {
      show <- c("id", "pop", "well", "het", "het.excess", "rare.burden",
                "rare.excess", "top.partner", "partner.well", "share",
                "kin.z", "flag")
      print(ind[ind$flag %in% c("suspect", "adjacent"), show],
            row.names = FALSE)
    }
  }
  if (verbose >= 3) {
    if (n.rare > 0) {
      cat(report("  Rare-only individuals (structure or mislabelling):\n"))
      print(ind[ind$flag == "rare-only",
                c("id", "pop", "het.excess", "rare.burden", "rare.excess")],
            row.names = FALSE)
    }
    if (!is.null(pairs)) {
      cat(report("  Plate-adjacent pairs with the largest residual kinship:\n"))
      print(utils::head(pairs, 10), row.names = FALSE)
    }
  }

  # PLOTS
  ind$class <- ifelse(ind$flag %in% c("suspect", "adjacent"),
                      "suspect", "other")
  lab <- ind[ind$class == "suspect", ]
  p1 <- ggplot2::ggplot(ind, ggplot2::aes(x = het, y = rare.burden,
                                          color = class)) +
    ggplot2::geom_point(size = 2, alpha = 0.8, na.rm = TRUE) +
    ggplot2::geom_text(data = lab, ggplot2::aes(label = id), size = 2.6,
                       vjust = -0.8, show.legend = FALSE, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(suspect = plot.colors[1],
                                           other = plot.colors[2]),
                                name = NULL) +
    ggplot2::labs(x = "Heterozygosity",
                  y = paste0("Rare-allele burden (freq < ", rare.freq,
                             " in own population)"),
                  title = "Contamination screen") +
    plot.theme
  if (plot.display) print(p1)
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p1, dir = plot.dir, file = plot.file,
                           verbose = verbose)
  }
  ind$class <- NULL

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN
  invisible(list(ind = ind, pairs = pairs, kinship = resid))
}
