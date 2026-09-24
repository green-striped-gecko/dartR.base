#' @name gl.report.contamination
#' @title Screens a genlight object for cross-contaminated samples
#' @family matched report
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
#' example "1-C2"; the well is the part after the last "-"), otherwise from \code{plate} and \code{well} columns.
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
#' Tier 3 asks whether a flagged individual's foreign alleles behave like a
#' minority of its reads or like half of them. A contaminant contributes a
#' small fraction of the DNA, so its alleles are called mainly at loci with
#' enough depth to catch a minority allele; a hybrid's alleles sit at 50
#' percent and are called at any depth. The test uses the loci fixed for one
#' allele in the individual's population (frequency below \code{rare.freq}
#' among its unflagged members) and for the other allele in its partner's
#' population (frequency at least 1 - \code{rare.freq} among its unflagged
#' members), so that every locus could show the foreign allele. Restricting
#' to such loci matters because deep loci are the conserved ones and would
#' otherwise dilute the rate; building the references from unflagged
#' animals matters because a second contaminated animal in the host
#' population would otherwise remove every locus at which it carries the
#' donor allele. \code{foreign.rate} is the proportion of these loci at
#' which the individual carries the foreign allele, the plain measure of
#' how much of the other population it holds.
#' These loci are split into quartiles of \code{rdepth}, the per-locus
#' average read depth that \code{gl.read.dart()} stores in
#' \code{@@other$loc.metrics}, and the proportion of them at which the
#' individual carries the foreign allele is compared across quartiles.
#' \code{depth.ratio} is that proportion in the deepest quartile over the
#' shallowest (with a half count added to each), \code{depth.p} a
#' chi-squared test for trend across the four quartiles, and \code{pattern}
#' is "saturated" when the foreign allele is already called at more than
#' half of the shallowest quartile (a heavy mixture or an F1 hybrid, which
#' depth cannot separate), "dose" when the ratio is at least
#' \code{depth.ratio} with p below 0.01 (a contaminant), "flat" when the
#' ratio is below \code{depth.flat} (a true heterozygote: hybrid, admixed
#' or mislabelled animal) and "unclear" between. The pattern is computed
#' for flagged individuals whose partner is in another population, with at
#' least 80 such loci and \code{min.share} foreign alleles among them, and
#' is NA otherwise or when the genlight has no \code{rdepth}; contamination
#' from the individual's own population, or from a population absent from
#' the dataset, cannot be tested this way. Read the pattern only when host
#' and partner are close enough for a hybrid to be plausible: across
#' species the plate-average depth is a poor proxy for the host's own depth,
#' and a "flat" call does not clear an animal whose foreign alleles are
#' fixed-absent in its species.
#' Locus average depth is a proxy for the individual's own read counts, so
#' the pattern does not change the flag; the decisive test remains allele
#' balance from per-individual read counts.
#'
#' Limits. The screen is a candidate list, not a verdict. Contamination
#' from the same population as the host raises only heterozygosity and
#' kinship, so a low rare-allele burden does not clear a suspect, and a
#' population whose members differ widely in heterozygosity (inbreeding,
#' admixture, relatives) can yield suspects that are not contaminated.
#' Power falls as host and contaminant become genetically closer: a mixture
#' between two populations of one species at a typical heterozygosity of
#' 0.3 shifts heterozygosity by less than the natural spread of many
#' populations. Tier 2 can only name a source that is in the dataset.
#'
#' The function densifies the genotype matrix with \code{as.matrix()}.
#'
#' The plot shows heterozygosity against rare-allele burden with suspects
#' labelled.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param rare.freq Allele frequency in the rest of the population below
#' which an allele counts as rare; must lie in (0, 0.5) [default 0.02].
#' @param min.n Minimum number of other individuals in the population called
#' at a locus for that locus to enter the rare-allele burden; at least 2
#' [default 5].
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
#' @param depth.ratio Foreign-allele rate in the deepest quartile of loci
#' over the shallowest at or above which, with trend p below 0.01, the
#' pattern is "dose" [default 1.5].
#' @param depth.flat The same ratio below which the pattern is "flat"
#' [default 1.25].
#' @param plate Data frame with columns id, plate and well (for example
#' "C4", in either case) giving plate positions; overrides positions found in the individual
#' metadata. Without it, positions come from \code{plate_location} in
#' \code{@@other$ind.metrics}, keyed by \code{service} when that column
#' exists, because a report that bundles orders repeats plate numbers
#' [default NULL].
#' @param plot.display If TRUE, resultant plots are displayed in the plot
#' window [default TRUE].
#' @param plot.theme A ggplot2 theme for the plot, for example
#' \code{theme_dartR()} [default theme_dartR()].
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
#' partner's sharing score), the tier 3 columns \code{foreign.rate},
#' \code{depth.ratio},
#' \code{depth.p} and \code{pattern}, and the flag;
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
#' @seealso \code{\link{gl.filter.contamination}},
#' \code{\link{gl.report.heterozygosity}},
#' \code{\link{gl.filter.heterozygosity}}
#' @importFrom stats mad median sd ave prop.trend.test
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
                                    depth.ratio = 1.5,
                                    depth.flat = 1.25,
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
    x <- .as_dartR(x)
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
  num.par <- list(rare.freq = rare.freq, min.n = min.n, z.flag = z.flag,
                  min.excess = min.excess, rare.min.excess = rare.min.excess,
                  min.share = min.share, share.tol = share.tol,
                  depth.ratio = depth.ratio, depth.flat = depth.flat)
  for (nm in names(num.par)) {
    v <- num.par[[nm]]
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v)) {
      stop(error(paste0("  ", nm, " must be a single finite number\n")))
    }
  }
  if (rare.freq <= 0 || rare.freq >= 0.5) {
    stop(error("  rare.freq must lie in (0, 0.5)\n"))
  }
  if (min.n < 2) {
    stop(error("  min.n must be at least 2\n"))
  }
  if (z.flag <= 0 || min.excess < 0 || rare.min.excess < 0) {
    stop(error("  z.flag must be positive and the excess margins ",
               "non-negative\n"))
  }
  if (min.share < 1 || share.tol < 0 || share.tol > 1) {
    stop(error("  min.share must be at least 1 and share.tol in [0, 1]\n"))
  }
  if (depth.flat < 1 || depth.ratio < depth.flat) {
    stop(error("  depth.flat must be at least 1 and depth.ratio at least ",
               "depth.flat\n"))
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
      # gl.read.dart() writes "<plate>-<well>"; split at the last "-" because
      # plate names may themselves contain "-"
      pl <- as.character(im$plate_location)
      dash <- grepl("-", pl)
      pos <- data.frame(plate = ifelse(dash, sub("-[^-]*$", "", pl), pl),
                        well = ifelse(dash, sub("^.*-", "", pl), NA))
    } else if (all(c("plate", "well") %in% names(im))) {
      pos <- data.frame(plate = im$plate, well = im$well)
    }
    # A report that bundles orders repeats plate numbers across orders
    if (!is.null(pos) && "service" %in% names(im)) {
      pos$plate <- paste(im$service, pos$plate, sep = "_")
    }
  }
  adjacent <- rep(NA, n.ind)
  pairs <- NULL
  w.ok <- FALSE
  if (!is.null(pos)) {
    pos$well <- toupper(trimws(as.character(pos$well)))
    w.row <- match(substr(pos$well, 1, 1), LETTERS)
    w.col <- suppressWarnings(as.integer(substring(pos$well, 2)))
    w.ok <- !is.na(w.row) & !is.na(w.col)
  }
  if (any(w.ok)) {
    is.adj <- function(i, j) {
      w.ok[i] && w.ok[j] &&
        identical(pos$plate[i], pos$plate[j]) &&
        (abs(w.row[i] - w.row[j]) + abs(w.col[i] - w.col[j])) == 1
    }
    adjacent <- vapply(seq_len(n.ind), function(i) is.adj(i, top.j[i]),
                       logical(1))
    # Adjacent pairs by direct lookup of the well below and the well to the
    # right of each individual, so each pair is found once without testing
    # all n(n - 1) / 2 pairs
    at <- split(which(w.ok), paste(pos$plate, w.row, w.col)[w.ok])
    nb <- lapply(which(w.ok), function(i) {
      c(at[[paste(pos$plate[i], w.row[i] + 1, w.col[i])]],
        at[[paste(pos$plate[i], w.row[i], w.col[i] + 1)]])
    })
    j <- unlist(nb, use.names = FALSE)
    i <- rep(which(w.ok), lengths(nb))
    pr <- cbind(pmin(i, j), pmax(i, j))
    pr <- pr[order(pr[, 2], pr[, 1]), , drop = FALSE]
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
    if (!is.null(pos) && any(!is.na(pos$well) & pos$well != "")) {
      cat(warn("  Warning: no plate well parses as a row letter and a",
               "column number (for example \"C4\"); adjacency not tested\n"))
    } else {
      cat(warn("  No plate positions found; adjacency not tested\n"))
    }
  }

  # Flags: heterozygosity is the required signal
  het.hit <- (het.z > z.flag & het.excess > min.excess) %in% TRUE
  rare.hit <- (rare.z > z.flag & rare.excess > rare.min.excess) %in% TRUE
  flag <- ifelse(het.hit, "suspect", ifelse(rare.hit, "rare-only", ""))
  flag[flag == "suspect" & adjacent %in% TRUE] <- "adjacent"

  # Tier 3: foreign-allele rate against locus depth, at loci fixed-different
  # between the individual's population and its partner's population
  depth <- x@other$loc.metrics$rdepth
  d.ratio <- rep(NA_real_, n.ind)
  d.p <- rep(NA_real_, n.ind)
  f.rate <- rep(NA_real_, n.ind)
  pattern <- rep(NA_character_, n.ind)
  if (is.null(depth) || !any(is.finite(depth))) {
    if (verbose >= 2) {
      cat(warn("  No rdepth in loc.metrics; depth pattern not tested\n"))
    }
  } else {
    if (verbose >= 2) {
      cat(report("  Tier 3: foreign-allele rate against locus depth\n"))
    }
    pop.freq <- function(idx) {
      n <- colSums(called[idx, , drop = FALSE])
      list(p = colSums(gm[idx, , drop = FALSE], na.rm = TRUE) / (2 * n), n = n)
    }
    # Reference sets are the unflagged members of each population: a second
    # contaminated animal in the host population would otherwise remove every
    # locus at which it carries the donor allele
    clean <- flag == ""
    for (i in which(!clean)) {
      j <- top.j[i]
      if (pops[j] == pops[i]) next
      own <- pop.freq(which(pops == pops[i] & clean))
      don <- pop.freq(which(pops == pops[j] & clean))
      g <- gm[i, ]
      ok <- called[i, ] & own$n >= min.n & don$n >= min.n & is.finite(depth)
      # inclusive on the donor side: an admixed or contaminated donor-population
      # animal carrying one host allele must not remove the locus
      alt.in <- ok & own$p < rare.freq & don$p >= 1 - rare.freq - 1e-9
      ref.in <- ok & own$p > 1 - rare.freq & don$p <= rare.freq + 1e-9
      cand <- which(alt.in | ref.in)
      if (length(cand) < 80) next
      foreign <- (alt.in & g >= 1)[cand] | (ref.in & g <= 1)[cand]
      f.rate[i] <- mean(foreign)
      if (sum(foreign) < min.share) next
      # rank-based quartiles so that tied depths cannot collapse a bin
      q <- ceiling(4 * rank(depth[cand], ties.method = "first") / length(cand))
      f <- tabulate(q[foreign], 4)
      n <- tabulate(q, 4)
      d.ratio[i] <- ((f[4] + 0.5) / (n[4] + 1)) / ((f[1] + 0.5) / (n[1] + 1))
      d.p[i] <- suppressWarnings(prop.trend.test(f, n)$p.value)
      # when the foreign allele is already called at most shallow loci, depth
      # no longer limits the calls and the ratio has no room to rise
      pattern[i] <- if (f[1] / n[1] > 0.5) {
        "saturated"
      } else if (d.ratio[i] >= depth.ratio && d.p[i] < 0.01) {
        "dose"
      } else if (d.ratio[i] < depth.flat) {
        "flat"
      } else {
        "unclear"
      }
    }
  }

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
                    foreign.rate = round(f.rate, 3),
                    depth.ratio = round(d.ratio, 2),
                    depth.p = signif(d.p, 2),
                    pattern = pattern,
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
                  "with an adjacent-well partner,",
                  sum(ind$pattern[ind$flag %in% c("suspect", "adjacent")]
                      %in% "dose"), "with a dose pattern,",
                  sum(ind$pattern[ind$flag %in% c("suspect", "adjacent")]
                      %in% "flat"), "flat );",
                  "rare-only:", n.rare, "\n"))
    if (n.sus > 0) {
      show <- c("id", "pop", "well", "het", "het.excess", "rare.burden",
                "rare.excess", "top.partner", "partner.well", "share",
                "kin.z", "foreign.rate", "depth.ratio", "pattern", "flag")
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
