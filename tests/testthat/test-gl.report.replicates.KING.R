# gl.report.replicates.KING: duplicate samples by KING-robust kinship and the
# heterozygote share, from genotype likelihoods, read counts or calls.

# --- Reference implementation (after Waples et al. 2019 and NgsRelate) -------

# Nine counts A-I of one pair from called genotypes; rows are individual 1
ref_calls9 <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  as.vector(table(factor(a[ok], 0:2), factor(b[ok], 0:2)))
}

# Nine expected counts of one pair by EM from linear-scale likelihoods
# (sites x 3), as in NgsRelate
ref_em9 <- function(L1, L2, tol = 1e-9, maxit = 5000) {
  g1 <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)
  g2 <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  M <- L1[, g1, drop = FALSE] * L2[, g2, drop = FALSE]
  p <- rep(1 / 9, 9)
  for (it in seq_len(maxit)) {
    pn <- p * colMeans(M / drop(M %*% p))
    done <- max(abs(pn - p)) < tol
    p <- pn
    if (done) break
  }
  p * nrow(M)
}

ref_stats <- function(v) {
  A <- v[1]; B <- v[2]; C <- v[3]; D <- v[4]; E <- v[5]
  F <- v[6]; G <- v[7]; H <- v[8]
  c(KING = (E - 2 * (C + G)) / (B + D + H + F + 2 * E),
    R0 = (C + G) / E,
    R1 = E / (B + C + D + F + G + H),
    het_share = E / min(B + E + H, D + E + F))
}

# Allele frequency by EM from likelihoods; L: individuals x sites x 3
ref_freq <- function(L, tol = 1e-8, maxit = 200) {
  f <- rep(0.25, dim(L)[2])
  for (it in seq_len(maxit)) {
    w0 <- t(t(L[, , 1]) * (1 - f)^2)
    w1 <- t(t(L[, , 2]) * 2 * f * (1 - f))
    w2 <- t(t(L[, , 3]) * f^2)
    fn <- colMeans((w1 + 2 * w2) / (w0 + w1 + w2), na.rm = TRUE) / 2
    done <- max(abs(fn - f)) < tol
    f <- fn
    if (done) break
  }
  f
}

# --- Fixtures ------------------------------------------------------------------

quiet <- function(expr) {
  capture.output(res <- suppressWarnings(expr))
  res
}

# platypus with one individual duplicated; 60% of the copy's heterozygotes
# are called homozygous (heterozygote dropout)
make_dropout <- function(frac = 0.6) {
  x <- platypus.gl
  G <- as.matrix(x)
  set.seed(1)
  g <- G["T27", ]
  h <- which(g == 1)
  k <- sample(h, round(frac * length(h)))
  g[k] <- sample(c(0, 2), length(k), replace = TRUE)
  dm <- rbind(G, g)
  rownames(dm) <- c(rownames(G), "dup")
  y <- new("genlight", gen = dm, ploidy = 2)
  indNames(y) <- rownames(dm)
  locNames(y) <- colnames(G)
  pop(y) <- factor(c(as.character(pop(x)), "TENTERFIELD"))
  suppressWarnings(gl.compliance.check(y, verbose = 0))
}

# A vcf of 8 samples: S1-S6 unrelated, S7 a deep duplicate of S1 and S8 a
# shallow duplicate of S2, with AD and PL from simulated reads
write_king_vcf <- function(n.loc = 400, seed = 1) {
  set.seed(seed)
  p <- runif(n.loc, 0.1, 0.9)
  geno <- t(replicate(6, rbinom(n.loc, 2, p)))
  geno <- rbind(geno, geno[1, ], geno[2, ])
  n <- nrow(geno)
  depth <- matrix(rpois(n * n.loc, 25), n)
  depth[8, ] <- rpois(n.loc, 2)
  e <- 0.01
  pa <- matrix(c(e, 0.5, 1 - e)[geno + 1], n)
  alt <- matrix(rbinom(length(depth), depth, pa), n)
  ref <- depth - alt
  pl <- lapply(c(e, 0.5, 1 - e), function(q) {
    round(-10 * (ref * log10(1 - q) + alt * log10(q)))
  })
  low <- do.call(pmin, pl)
  pl <- lapply(pl, function(m) m - low)
  gt <- matrix(c("0/0", "0/1", "1/1")[max.col(-cbind(
    as.vector(pl[[1]]), as.vector(pl[[2]]), as.vector(pl[[3]])),
    ties.method = "first")], n)
  gt[depth == 0] <- "./."
  cell <- matrix(paste0(gt, ":", depth, ":", ref, ",", alt, ":", pl[[1]], ",",
                        pl[[2]], ",", pl[[3]]), n)
  recs <- vapply(seq_len(n.loc), function(j) {
    paste(c("chr1", j * 10, paste0("snp", j), "A", "G", "50", "PASS", ".",
            "GT:DP:AD:PL", cell[, j]), collapse = "\t")
  }, character(1))
  f <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=100000>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", paste0("S", 1:n)), collapse = "\t"),
    recs), f)
  f
}

read_king_vcf <- function() {
  skip_if_not_installed("vcfR")
  quiet(gl.read.vcf(write_king_vcf(), read.counts = TRUE, likelihoods = TRUE,
                    verbose = 0))
}

pair_row <- function(tab, a, b) {
  tab[(tab$ind1 == a & tab$ind2 == b) | (tab$ind1 == b & tab$ind2 == a), ]
}

# --- Calls ---------------------------------------------------------------------

test_that("calls: pair statistics equal a direct computation", {
  x <- make_dropout()
  res <- quiet(gl.report.replicates.KING(x, plot.out = FALSE, verbose = 0))
  expect_identical(res$input, "calls")
  G <- as.matrix(x)
  f <- colMeans(G, na.rm = TRUE) / 2
  use <- which(f >= 0.05 & f <= 0.95)
  expect_equal(res$parameters$n.sites, length(use))
  for (pr in list(c("T27", "dup"), c("T27", "T35"), c("SDS4", "SUS20"))) {
    v <- ref_calls9(G[pr[1], use], G[pr[2], use])
    row <- pair_row(res$table.all, pr[1], pr[2])
    expect_equal(nrow(row), 1)
    expect_equal(row$nloc, sum(v))
    expect_equal(unlist(row[, c("KING", "R0", "R1", "het_share")]),
                 ref_stats(v), ignore_attr = TRUE)
  }
  # every pair once, ind1 earlier in the object than ind2
  n <- nInd(x)
  expect_equal(nrow(res$table.all), n * (n - 1) / 2)
  expect_true(all(match(res$table.all$ind1, indNames(x)) <
                    match(res$table.all$ind2, indNames(x))))
})

test_that("calls: a duplicate with heterozygote dropout is flagged by het_share", {
  x <- make_dropout()
  # the identity screen misses it
  old <- quiet(gl.report.replicates(x, plot.out = FALSE, verbose = 0))
  expect_false(any(pair_row(old$table.rep, "T27", "dup")$perc > 0.95))

  res <- quiet(gl.report.replicates.KING(x, plot.out = FALSE, verbose = 0))
  row <- pair_row(res$table.all, "T27", "dup")
  expect_lt(row$KING, 0.354)
  expect_gt(row$het_share, 0.95)
  expect_true(row$flagged)
  # flagged only by het_share, but classed duplicate
  expect_identical(row$class, "duplicate")
  expect_true(nrow(pair_row(res$table.rep, "T27", "dup")) == 1)
  # the copy has the same missing genotypes as T27, so the name decides
  keep <- c("T27", "dup")[order(res$ind.miss[c("T27", "dup")], c("T27", "dup"),
                                method = "radix")[1]]
  # platypus has two duplicate pairs of its own (T3-T5, T28-T42)
  expect_true(setdiff(c("T27", "dup"), keep) %in% res$ind.list.drop)
  expect_false(keep %in% res$ind.list.drop)
  expect_setequal(res$ind.list.drop, c("T3", "T42", setdiff(c("T27", "dup"), keep)))

  # without the heterozygote share the pair is not flagged
  res2 <- quiet(gl.report.replicates.KING(x, het_threshold = 1.01,
                                          plot.out = FALSE, verbose = 0))
  expect_false(pair_row(res2$table.all, "T27", "dup")$flagged)
})

test_that("calls: the report warns that KING underestimates duplicates", {
  x <- make_dropout()
  out <- capture.output(suppressWarnings(
    gl.report.replicates.KING(x, plot.out = FALSE, verbose = 1)))
  expect_true(any(grepl("heterozygote dropout", out)))
})

test_that("KING classes follow the Manichaikul bins", {
  x <- make_dropout()
  res <- quiet(gl.report.replicates.KING(x, plot.out = FALSE, verbose = 0))
  tab <- res$table.all
  ok <- !tab$flagged & !is.na(tab$KING)
  expected <- cut(tab$KING[ok], c(-Inf, 0.0442, 0.0884, 0.177, Inf),
                  labels = c("unrelated", "third degree", "second degree",
                             "first degree"))
  expect_identical(tab$class[ok], as.character(expected))
  expect_true(all(tab$class[tab$flagged] == "duplicate"))
})

test_that("pairs with fewer than loc_threshold shared sites are never flagged", {
  x <- make_dropout()
  res <- quiet(gl.report.replicates.KING(x, loc_threshold = 5000,
                                         plot.out = FALSE, verbose = 0))
  expect_true(all(is.na(res$table.all$KING)))
  expect_true(all(is.na(res$table.all$class)))
  expect_false(any(res$table.all$flagged))
  expect_equal(nrow(res$table.rep), 0)
  expect_length(res$ind.list.drop, 0)
  expect_length(res$ind.list.rep, 0)
  expect_true(all(res$table.all$nloc > 0))
})

test_that("the report returns the documented elements", {
  x <- make_dropout()
  res <- quiet(gl.report.replicates.KING(x, plot.out = FALSE, verbose = 0))
  expect_named(res, c("table.rep", "ind.list.drop", "ind.list.rep",
                      "table.all", "ind.miss", "input", "parameters", "plot"),
               ignore.order = TRUE)
  expect_true(all(c("ind1", "ind2", "nloc", "KING", "R0", "R1", "het_share",
                    "class", "group", "ind1_miss", "ind2_miss",
                    "ind_to_drop") %in% names(res$table.rep)))
  expect_identical(names(res$ind.miss), indNames(x))
  expect_equal(unname(res$ind.miss), unname(rowMeans(is.na(as.matrix(x)))))
  expect_true(all(c("T27", "dup") %in% names(res$ind.list.rep)))
  expect_null(res$plot)
  capture.output(v <- withVisible(gl.report.replicates.KING(
    x, plot.out = FALSE, verbose = 0)))
  expect_false(v$visible)
})

test_that("the report works on FBM-backed objects", {
  x <- make_dropout()
  a <- quiet(gl.report.replicates.KING(x, plot.out = FALSE, verbose = 0))
  b <- quiet(gl.report.replicates.KING(gl.gen2fbm(x, verbose = 0),
                                       plot.out = FALSE, verbose = 0))
  expect_equal(a$table.all, b$table.all)
})

# --- Groups ----------------------------------------------------------------------

# B carries the alleles of both A and C, as a mixture of the two would: it
# shares every heterozygote of each, so both pairs are flagged, while A and
# C are unrelated
make_chain <- function() {
  x <- platypus.gl
  G <- as.matrix(x)
  a <- G["T27", ]
  c2 <- G["SUS20", ]
  has.ref <- (a < 2 | c2 < 2)
  has.alt <- (a > 0 | c2 > 0)
  b <- ifelse(has.ref & has.alt, 1, ifelse(has.alt, 2, 0))
  dm <- rbind(G, b)
  rownames(dm) <- c(rownames(G), "mix")
  y <- new("genlight", gen = dm, ploidy = 2)
  indNames(y) <- rownames(dm)
  locNames(y) <- colnames(G)
  pop(y) <- factor(c(as.character(pop(x)), "mix"))
  suppressWarnings(gl.compliance.check(y, verbose = 0))
}

test_that("a chain through a mixed sample is grouped and listed", {
  x <- make_chain()
  cc <- "SUS20"
  out <- capture.output(res <- suppressWarnings(
    gl.report.replicates.KING(x, plot.out = FALSE, verbose = 1)))
  expect_true(pair_row(res$table.all, "T27", "mix")$flagged)
  expect_true(pair_row(res$table.all, cc, "mix")$flagged)
  expect_false(pair_row(res$table.all, "T27", cc)$flagged)
  g <- res$table.rep$group[res$table.rep$ind1 %in% c("T27", cc, "mix") &
                             res$table.rep$ind2 %in% c("T27", cc, "mix")]
  expect_length(unique(g), 1)
  # the warning names the unflagged pair inside the group
  w <- out[grepl("not flagged", out) | grepl(cc, out, fixed = TRUE)]
  expect_true(any(grepl("not flagged", out)))
  expect_true(any(grepl("T27", w) & grepl(cc, w, fixed = TRUE)))
  # one sample of the three is kept: the one with the fewest missing
  # genotypes
  three <- c("T27", cc, "mix")
  keep <- three[order(res$ind.miss[three], three, method = "radix")[1]]
  expect_setequal(intersect(res$ind.list.drop, three), setdiff(three, keep))
})

test_that("no flagged pairs returns empty results", {
  x <- platypus.gl[1:20, ]
  res <- quiet(gl.report.replicates.KING(x, king_threshold = 2,
                                         het_threshold = 2,
                                         plot.out = FALSE, verbose = 0))
  expect_equal(nrow(res$table.rep), 0)
  expect_length(res$ind.list.drop, 0)
  expect_length(res$ind.list.rep, 0)
  expect_equal(nrow(res$table.all), 190)
})

# --- Likelihoods and read counts -----------------------------------------------

test_that("likelihoods: statistics equal the NgsRelate EM for every pair", {
  skip_if_not_installed("Rcpp")
  skip_if_not_installed("RcppParallel")
  x <- read_king_vcf()
  res <- quiet(gl.report.replicates.KING(x, loc_threshold = 50,
                                         plot.out = FALSE, verbose = 0))
  expect_identical(res$input, "likelihoods")
  gm <- x@other$geno.metrics
  L <- array(c(10^(-gm$pl.ref / 10), 10^(-gm$pl.het / 10),
               10^(-gm$pl.alt / 10)), c(nInd(x), nLoc(x), 3))
  f <- ref_freq(L)
  use <- which(f >= 0.05 & f <= 0.95)
  expect_equal(res$parameters$n.sites, length(use))
  L <- L[, use, , drop = FALSE]
  has <- !is.na(L[, , 1])
  tab <- res$table.all
  for (r in seq_len(nrow(tab))) {
    i <- match(tab$ind1[r], indNames(x))
    j <- match(tab$ind2[r], indNames(x))
    s <- has[i, ] & has[j, ]
    v <- ref_em9(matrix(L[i, s, ], ncol = 3), matrix(L[j, s, ], ncol = 3))
    expect_equal(tab$nloc[r], sum(s))
    expect_equal(unlist(tab[r, c("KING", "R0", "R1", "het_share")]),
                 ref_stats(v), tolerance = 1e-10, ignore_attr = TRUE)
  }
  # both duplicates are flagged, and nothing else
  fl <- tab[tab$flagged, ]
  expect_setequal(paste(fl$ind1, fl$ind2), c("S1 S7", "S2 S8"))
})

test_that("counts: statistics equal the EM on binomial read likelihoods", {
  skip_if_not_installed("Rcpp")
  skip_if_not_installed("RcppParallel")
  x <- read_king_vcf()
  e <- 0.001
  res <- quiet(gl.report.replicates.KING(x, input = "counts", error_rate = e,
                                         loc_threshold = 50,
                                         plot.out = FALSE, verbose = 0))
  expect_identical(res$input, "counts")
  G <- as.matrix(x)
  f <- colMeans(G, na.rm = TRUE) / 2
  use <- which(f >= 0.05 & f <= 0.95)
  R <- x@other$geno.metrics$ref.count[, use]
  A <- x@other$geno.metrics$alt.count[, use]
  L <- array(c((1 - e)^R * e^A, 0.5^(R + A), e^R * (1 - e)^A),
             c(nInd(x), length(use), 3))
  L[rep(R + A == 0, 3)] <- NA
  has <- !is.na(L[, , 1])
  tab <- res$table.all
  for (r in seq_len(nrow(tab))) {
    i <- match(tab$ind1[r], indNames(x))
    j <- match(tab$ind2[r], indNames(x))
    s <- has[i, ] & has[j, ]
    v <- ref_em9(matrix(L[i, s, ], ncol = 3), matrix(L[j, s, ], ncol = 3))
    expect_equal(tab$nloc[r], sum(s))
    expect_equal(unlist(tab[r, c("KING", "R0", "R1", "het_share")]),
                 ref_stats(v), tolerance = 1e-8, ignore_attr = TRUE)
  }
  fl <- tab[tab$flagged, ]
  expect_setequal(paste(fl$ind1, fl$ind2), c("S1 S7", "S2 S8"))
})

test_that("input auto picks likelihoods, then counts, then calls", {
  skip_if_not_installed("Rcpp")
  skip_if_not_installed("RcppParallel")
  x <- read_king_vcf()
  out <- capture.output(res <- suppressWarnings(
    gl.report.replicates.KING(x, loc_threshold = 50, plot.out = FALSE,
                              verbose = 2)))
  expect_identical(res$input, "likelihoods")
  expect_true(any(grepl("likelihoods", out)))
  y <- x
  y@other$geno.metrics[c("pl.ref", "pl.het", "pl.alt")] <- NULL
  expect_identical(quiet(gl.report.replicates.KING(
    y, loc_threshold = 50, plot.out = FALSE, verbose = 0))$input, "counts")
  y@other$geno.metrics <- NULL
  expect_identical(quiet(gl.report.replicates.KING(
    y, loc_threshold = 50, plot.out = FALSE, verbose = 0))$input, "calls")
  # asking for an input the object does not hold is an error
  expect_error(quiet(gl.report.replicates.KING(y, input = "likelihoods",
                                               verbose = 0)), "likelihoods")
  expect_error(quiet(gl.report.replicates.KING(y, input = "counts",
                                               verbose = 0)), "counts")
})

test_that("a shallow duplicate scores higher from likelihoods than from calls", {
  skip_if_not_installed("Rcpp")
  skip_if_not_installed("RcppParallel")
  x <- read_king_vcf()
  pl <- quiet(gl.report.replicates.KING(x, loc_threshold = 50,
                                        plot.out = FALSE, verbose = 0))
  ca <- quiet(gl.report.replicates.KING(x, input = "calls", loc_threshold = 50,
                                        plot.out = FALSE, verbose = 0))
  expect_gt(pair_row(pl$table.all, "S2", "S8")$KING,
            pair_row(ca$table.all, "S2", "S8")$KING)
})

# --- Plot ------------------------------------------------------------------------

test_that("the plot is returned and printed only when asked", {
  x <- make_dropout()
  pdf(NULL)
  on.exit(dev.off())
  res <- quiet(gl.report.replicates.KING(x, verbose = 1))
  expect_s3_class(res$plot, "patchwork")
  expect_s3_class(res$plot, "ggplot")
  # at verbose 0 nothing is printed, but the plot is still returned
  t0 <- tempfile(fileext = ".pdf")
  pdf(t0)
  dev.off()
  tf <- tempfile(fileext = ".pdf")
  pdf(tf)
  res0 <- quiet(gl.report.replicates.KING(x, verbose = 0))
  dev.off()
  expect_lt(file.info(tf)$size, file.info(t0)$size + 500)
  expect_s3_class(res0$plot, "patchwork")
})
