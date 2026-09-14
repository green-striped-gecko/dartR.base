# Characterization tests for gl.impute
# Baseline snapshotted before review (review-gl.impute, dev at ddaed27).
# Phase C (2026-09-06): assertions formerly marked [baseline defect Fn] have
# been flipped to the approved behaviour and are marked [approved Fn].

# Small synthetic genlight: rows = individuals, 0/1/2/NA
mkgl <- function(mat, pops) {
  g <- methods::new("genlight", mat, ploidy = 2)
  adegenet::indNames(g) <- paste0("i", seq_len(nrow(mat)))
  adegenet::locNames(g) <- paste0("L", seq_len(ncol(mat)))
  adegenet::pop(g) <- factor(pops)
  g
}

# Small synthetic presence-absence genlight: 0/1/NA, ploidy 1
mkgs <- function(mat, pops) {
  g <- methods::new("genlight", mat, ploidy = 1)
  adegenet::indNames(g) <- paste0("i", seq_len(nrow(mat)))
  adegenet::locNames(g) <- paste0("L", seq_len(ncol(mat)))
  adegenet::pop(g) <- factor(pops)
  g
}

test_that("frequency and HW impute per population (deterministic at q=0/1)", {
  pdf(NULL); on.exit(dev.off())
  m <- rbind(matrix(0, 6, 10), matrix(2, 6, 10))
  m[, 10] <- rep(c(0, 2), 6)   # keep one polymorphic locus
  na_cells <- cbind(c(1, 2, 3, 7, 8, 9), c(1, 3, 5, 2, 4, 6))
  m[na_cells] <- NA
  g <- mkgl(m, rep(c("A", "B"), each = 6))
  for (meth in c("frequency", "HW")) {
    set.seed(42)
    out <- as.matrix(gl.impute(g, method = meth, verbose = 0))
    # popA fixed for 0 (q=0), popB fixed for 2 (q=1): imputation is per-pop
    expect_equal(unname(out[na_cells]), c(0, 0, 0, 2, 2, 2))
    # originally scored cells untouched
    expect_true(all(out[!is.na(m)] == m[!is.na(m)]))
    expect_equal(dim(out), dim(m))
    expect_equal(sum(is.na(out)), 0)
  }
})

test_that("frequency is a deterministic expected-dosage fill, distinct from HW", {
  # [approved F1] 'frequency' now fills each missing genotype with the
  # expected dosage 2q (per-population) rounded to the nearest valid
  # genotype, ties (2q = 0.5 or 1.5) to the heterozygote. Previously it
  # drew two alleles Bernoulli(q), distributionally identical to 'HW'.
  pdf(NULL); on.exit(dev.off())
  set.seed(7)
  m <- matrix(rep(c(0, 2), each = 10), nrow = 20, ncol = 200)
  miss <- matrix(runif(length(m)) < 0.3, nrow = 20)
  m[miss] <- NA
  g <- mkgl(m, rep("A", 20))
  set.seed(99)
  out <- as.matrix(gl.impute(g, method = "frequency", verbose = 0))
  # hand-computed statistic, every imputed cell: q from the scored cells
  qcol <- colMeans(m, na.rm = TRUE) / 2
  d <- 2 * qcol
  expected <- ifelse(d < 0.5, 0, ifelse(d <= 1.5, 1, 2))
  exp_mat <- matrix(expected, nrow = nrow(m), ncol = ncol(m), byrow = TRUE)
  expect_equal(unname(out[miss]), unname(exp_mat[miss]))  # [approved F1]
  # deterministic: same seed and different seed give identical results
  set.seed(99)
  a <- as.matrix(gl.impute(g, method = "frequency", verbose = 0))
  set.seed(12345)
  b <- as.matrix(gl.impute(g, method = "frequency", verbose = 0))
  expect_identical(a, b)                                  # [approved F1]
  # HW on the same fixture still scatters 0/1/2 -> the methods now differ
  set.seed(99)
  hw <- as.matrix(gl.impute(g, method = "HW", verbose = 0))[miss]
  expect_true(all(c(0, 1, 2) %in% hw))
})

test_that("frequency fill is per-population, with the documented tie rule", {
  # [approved F1] cell-level check on a two-population fixture with known
  # per-population allele frequencies, including both tie cases.
  pdf(NULL); on.exit(dev.off())
  A <- rbind(c(0, 2, 1, 0, 2, 1, 2, 0),   # masked row (values irrelevant)
             c(0, 2, 1, 0, 2, 1, 2, 0),
             c(0, 2, 1, 0, 2, 1, 2, 1),
             c(0, 2, 1, 0, 2, 0, 1, 1),
             c(0, 2, 1, 1, 1, 0, 1, 2))
  B <- rbind(c(2, 0, 0, 2, 0, 2, 0, 2),   # masked row (values irrelevant)
             c(2, 0, 0, 2, 0, 2, 0, 2),
             c(2, 0, 0, 2, 0, 2, 0, 2),
             c(2, 0, 0, 2, 1, 2, 0, 2),
             c(2, 0, 1, 2, 1, 0, 0, 2))
  m <- rbind(A, B)
  m[1, ] <- NA                            # individual 1 (popA) all missing
  m[6, ] <- NA                            # individual 6 (popB) all missing
  g <- mkgl(m, rep(c("A", "B"), each = 5))
  set.seed(11)
  out <- as.matrix(gl.impute(g, method = "frequency", verbose = 0))
  # popA scored cells per locus give 2q = 0, 2, 1, 0.25, 1.75, 0.5, 1.5, 1
  expect_equal(unname(out[1, ]), c(0, 2, 1, 0, 2, 1, 1, 1))  # [approved F1]
  # popB scored cells per locus give 2q = 2, 0, 0.25, 2, 0.5, 1.5, 0, 2
  expect_equal(unname(out[6, ]), c(2, 0, 0, 2, 1, 1, 0, 2))  # [approved F1]
  # scored cells untouched
  expect_true(all(out[!is.na(m)] == m[!is.na(m)]))
})

test_that("neighbour copies from the nearest neighbour, then the next", {
  pdf(NULL); on.exit(dev.off())
  m <- rbind(c(NA, rep(1, 9)),
             c(2,  rep(1, 9)),      # identical elsewhere -> nearest
             c(0,  rep(0, 9)),
             c(0,  rep(2, 9)),
             c(0,  rep(2, 9)))
  g <- mkgl(m, rep("A", 5))
  out <- as.matrix(gl.impute(g, method = "neighbour", verbose = 0))
  expect_equal(unname(out[1, 1]), 2)             # from i2, not the distant 0s
  # nearest also NA at the target locus -> falls through to second nearest
  m2 <- m
  m2[2, 1] <- NA
  m2[3, ] <- c(0, rep(1, 8), 0)                  # i3 now second nearest to i1
  g2 <- mkgl(m2, rep("A", 5))
  out2 <- as.matrix(gl.impute(g2, method = "neighbour",
                              fill.residual = FALSE, verbose = 0))
  expect_equal(unname(out2[1, 1]), 0)
  # deterministic across runs
  r1 <- as.matrix(gl.impute(g, method = "neighbour", verbose = 0))
  r2 <- as.matrix(gl.impute(g, method = "neighbour", verbose = 0))
  expect_identical(r1, r2)
})

test_that("random draws uniformly from 0:2 and is seed-reproducible", {
  pdf(NULL); on.exit(dev.off())
  set.seed(7)
  m <- matrix(rep(c(0, 2), each = 10), nrow = 20, ncol = 200)
  miss <- matrix(runif(length(m)) < 0.3, nrow = 20)
  m[miss] <- NA
  g <- mkgl(m, rep("A", 20))
  set.seed(5)
  imp <- as.matrix(gl.impute(g, method = "random", verbose = 0))[miss]
  expect_true(all(imp %in% 0:2))
  expect_gt(mean(imp == 1), 0.25)                # uniform ~1/3, not HW 1/2
  expect_lt(mean(imp == 1), 0.42)
  set.seed(77)
  a <- as.matrix(gl.impute(g, method = "random", verbose = 0))
  set.seed(77)
  b <- as.matrix(gl.impute(g, method = "random", verbose = 0))
  expect_identical(a, b)
})

test_that("residual-NA contract: pop-all-NA vs globally-all-NA loci", {
  pdf(NULL); on.exit(dev.off())
  m <- rbind(matrix(0, 4, 5), matrix(2, 4, 5))
  m[, 5] <- rep(c(0, 2), 4)                      # polymorphic guard locus
  m[1:4, 3] <- NA                                # L3 all-NA in popA only
  g <- mkgl(m, rep(c("A", "B"), each = 4))
  oF <- gl.impute(g, method = "frequency", fill.residual = FALSE, verbose = 0)
  expect_true(all(is.na(as.matrix(oF)[1:4, 3])))
  set.seed(3)
  oT <- gl.impute(g, method = "frequency", fill.residual = TRUE, verbose = 0)
  # filled from the GLOBAL profile (q = 1 across popB) -> all 2s
  expect_equal(unname(as.matrix(oT)[1:4, 3]), rep(2, 4))
  # a locus NA across ALL individuals stays NA even with fill.residual = TRUE
  m2 <- m
  m2[, 4] <- NA
  g2 <- mkgl(m2, rep(c("A", "B"), each = 4))
  set.seed(3)
  oG <- gl.impute(g2, method = "frequency", verbose = 0)
  expect_true(all(is.na(as.matrix(oG)[, 4])))
})

test_that("testset.gl: metadata resynced, history appended, silence at verbose 0", {
  pdf(NULL); on.exit(dev.off())
  gl <- gl.compliance.check(testset.gl[1:40, 1:80], verbose = 0)
  m0 <- as.matrix(gl)
  allna <- which(colSums(!is.na(m0)) == 0)
  expect_length(allna, 2)                        # dartR.data 1.2.5
  set.seed(1)
  o <- capture.output(out <- gl.impute(gl, method = "neighbour", verbose = 0))
  # Known sibling leak at ddaed27: utils.dist.ind.snp prints its banner
  # ungated at any verbosity (fix pending on review-utils.dist.ind.snp);
  # gl.impute itself adds nothing at verbose 0.
  o <- o[!grepl("Calculating the .*distance matrix", o)]
  expect_length(o, 0)                            # VRB5: fully silent
  m1 <- as.matrix(out)
  expect_identical(dim(m1), dim(m0))
  expect_identical(indNames(out), indNames(gl))
  expect_identical(locNames(out), locNames(gl))
  expect_true(all(m1[!is.na(m0)] == m0[!is.na(m0)]))
  # loc.metrics recalculated by the internal compliance pass (DAT4)
  cr <- out@other$loc.metrics$CallRate
  expect_true(all(cr[-allna] == 1))
  expect_equal(unname(cr[allna]), c(0, 0))
  # residual NAs are exactly the all-NA loci
  expect_equal(sum(is.na(m1)), length(allna) * nInd(gl))
  # per-individual metadata still aligned (DAT2)
  expect_identical(as.character(out@other$ind.metrics$id), indNames(out))
  expect_identical(as.character(pop(out)), as.character(pop(gl)))
  expect_true(all(ploidy(out) == 2))
  # history: input's entries carried over, gl.impute appended (FS8)
  h <- out@other$history
  expect_equal(length(h), length(gl@other$history) + 1)
  expect_equal(deparse(h[[length(h)]][[1]]), "gl.impute")
})

test_that("frequency path preserves individual order across seppop/rbind", {
  pdf(NULL); on.exit(dev.off())
  gl <- gl.compliance.check(testset.gl[1:40, 1:80], verbose = 0)
  m0 <- as.matrix(gl)
  set.seed(2)
  out <- gl.impute(gl, method = "frequency", verbose = 0)
  expect_identical(indNames(out), indNames(gl))
  m1 <- as.matrix(out)
  expect_true(all(m1[!is.na(m0)] == m0[!is.na(m0)]))
})

test_that("no missing data: genotypes pass through unchanged", {
  pdf(NULL); on.exit(dev.off())
  set.seed(8)
  m <- matrix(sample(0:2, 40, TRUE), 4, 10)
  g <- mkgl(m, rep(c("A", "B"), each = 2))
  out <- gl.impute(g, method = "frequency", verbose = 0)
  expect_identical(unname(as.matrix(out)), unname(m))
})

test_that("invalid method fails with an informative error", {
  # [approved F4] method is now validated before any work; a typo
  # previously fell through every branch and died with the internal
  # error "object 'x3' not found".
  pdf(NULL); on.exit(dev.off())
  set.seed(8)
  m <- matrix(sample(0:2, 40, TRUE), 4, 10)
  g <- mkgl(m, rep(c("A", "B"), each = 2))
  expect_error(gl.impute(g, method = "frequencies", verbose = 0),
               "method must be one of")          # [approved F4]
})

test_that("all-missing warning fires only for genuinely all-missing loci", {
  # [approved F5] detection is now glNA(y) == ploidy * nInd(y); previously
  # glNA(y) > nInd(y) reported any locus with >50% missing as all-missing.
  pdf(NULL); on.exit(dev.off())
  m <- matrix(rep(c(0, 2), 12), 6, 4)
  m[1:4, 2] <- NA                                # 67% missing, NOT all
  g <- mkgl(m, rep("A", 6))
  o <- capture.output(out <- gl.impute(g, method = "frequency", verbose = 2))
  expect_false(any(grepl("loci with all missing values", o)))  # [approved F5]
  # a genuinely all-missing locus is still reported, with the right count
  m2 <- m
  m2[, 2] <- NA
  g2 <- mkgl(m2, rep("A", 6))
  o2 <- capture.output(out2 <- gl.impute(g2, method = "frequency", verbose = 3))
  expect_true(any(grepl("has\\s+1\\s+loci with all missing values", o2)))
  expect_true(any(grepl("0 values to be imputed", o2)))
})

test_that("per-population warnings name every affected population", {
  # [approved F6] the random branch printed its warning after the
  # population loop, so only the last population's numbers were ever
  # reported; a pop with all-missing loci earlier in the list was silent.
  pdf(NULL); on.exit(dev.off())
  m <- matrix(rep(c(0, 2), 16), 8, 4)
  m[1:4, 2] <- NA                                # L2 all-NA in popA only
  g <- mkgl(m, rep(c("A", "B"), each = 4))
  set.seed(9)
  o <- capture.output(out <- gl.impute(g, method = "random", verbose = 2))
  hit <- grep("loci with all missing values", o, value = TRUE)
  expect_length(hit, 1)
  expect_true(grepl("Population\\s+A\\s", hit))  # [approved F6]
})

test_that("SilicoDArT: every method keeps the 0/1 domain", {
  # [approved F3] the samplers are now domain-aware: 'random' draws from
  # 0:1, 'frequency' is a deterministic majority-band fill, 'HW' is a
  # Bernoulli draw with the band frequency. Previously 'random' and
  # 'frequency' wrote genotype 2s into 0/1 data and 'frequency'/'HW'
  # crashed in the residual pass with "negative probability".
  pdf(NULL); on.exit(dev.off())
  gs <- gl.compliance.check(testset.gs[1:30, 1:60], verbose = 0)
  m0 <- as.matrix(gs)
  expect_true(all(m0 %in% c(0, 1, NA)))
  for (meth in c("random", "frequency", "HW", "neighbour")) {
    set.seed(4)
    r <- gl.impute(gs, method = meth, verbose = 0)
    mr <- as.matrix(r)
    expect_true(all(mr %in% c(0, 1, NA)), label = meth)   # [approved F3]
    expect_true(all(mr[!is.na(m0)] == m0[!is.na(m0)]), label = meth)
    expect_true(all(ploidy(r) == 1), label = meth)
  }
  # 'frequency' imputation is deterministic across seeds (the optional
  # residual fill remains a documented random draw, so it is off here)
  set.seed(21)
  f1 <- as.matrix(gl.impute(gs, method = "frequency",
                            fill.residual = FALSE, verbose = 0))
  set.seed(22)
  f2 <- as.matrix(gl.impute(gs, method = "frequency",
                            fill.residual = FALSE, verbose = 0))
  expect_identical(f1, f2)                       # [approved F1/F3]
  # 'HW' and 'random' are seed-reproducible
  for (meth in c("HW", "random")) {
    set.seed(23)
    a <- as.matrix(gl.impute(gs, method = meth, verbose = 0))
    set.seed(23)
    b <- as.matrix(gl.impute(gs, method = meth, verbose = 0))
    expect_identical(a, b, label = meth)
  }
  # verbose 0 remains fully silent on presence-absence data
  set.seed(4)
  o <- capture.output(rs <- gl.impute(gs, method = "frequency", verbose = 0))
  expect_length(o, 0)
})

test_that("SilicoDArT frequency fill is the per-population majority band", {
  # [approved F3] cell-level check with known band frequencies, including
  # the tie (band frequency exactly 0.5 resolves to presence).
  pdf(NULL); on.exit(dev.off())
  m <- rbind(c(NA, NA, NA),
             c(1, 0, 1),
             c(1, 0, 0),
             c(1, 0, 0),
             c(1, 1, 0))
  # scored band frequencies per locus: 1.0, 0.25, 0.25 -> fills 1, 0, 0
  g <- mkgs(m, rep("A", 5))
  set.seed(31)
  out <- as.matrix(gl.impute(g, method = "frequency", verbose = 0))
  expect_equal(unname(out[1, ]), c(1, 0, 0))     # [approved F3]
  # tie: two of four scored bands present -> q = 0.5 -> presence (1)
  m2 <- rbind(c(NA, NA), c(1, 1), c(1, 0), c(0, 1), c(0, 0))
  g2 <- mkgs(m2, rep("A", 5))
  out2 <- as.matrix(gl.impute(g2, method = "frequency", verbose = 0))
  expect_equal(unname(out2[1, ]), c(1, 1))       # [approved F3] tie -> 1
})

test_that("beagle: informative guards without binaries", {
  # [approved F3, F7] beagle is blocked for presence-absence data, and the
  # jar/java/PLINK toolchain is checked up front with stop(error()) rather
  # than surfacing raw system errors (or returning -1 for missing R.utils).
  pdf(NULL); on.exit(dev.off())
  gs <- gl.compliance.check(testset.gs[1:10, 1:20], verbose = 0)
  expect_error(gl.impute(gs, method = "beagle", verbose = 0),
               "SNP genotype data only")         # [approved F3]
  skip_if_not_installed("R.utils")
  set.seed(8)
  m <- matrix(sample(0:2, 40, TRUE), 4, 10)
  g <- mkgl(m, rep(c("A", "B"), each = 2))
  nojar <- file.path(tempdir(), "gl-impute-empty-dir")
  dir.create(nojar, showWarnings = FALSE)
  expect_error(gl.impute(g, method = "beagle",
                         beagle.bin.path = nojar, verbose = 0),
               "cannot find the beagle jar")     # [approved F7]
})

test_that("FBM-backed run is identical to the dense run, residual NAs preserved", {
  # [approved F2] residual NAs are now written to the FBM backing as
  # code 3 (NA under CODE_012, as in gl.gen2fbm); previously the raw
  # assignment coerced NA to 0, fabricating homozygous-reference
  # genotypes at the all-NA loci.
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off())
  gl <- gl.compliance.check(testset.gl[1:40, 1:80], verbose = 0)
  fb <- suppressWarnings(gl.gen2fbm(gl, verbose = 0))
  expect_equal(sum(is.na(as.matrix(fb))), sum(is.na(as.matrix(gl))))
  set.seed(5)
  rfb <- suppressWarnings(gl.impute(fb, method = "neighbour", verbose = 0))
  set.seed(5)
  rd <- gl.impute(gl, method = "neighbour", verbose = 0)
  md <- as.matrix(rd)
  mr <- as.matrix(rfb)
  na_dense <- which(is.na(md))
  expect_gt(length(na_dense), 0)
  expect_true(all(is.na(mr[na_dense])))          # [approved F2]
  # cell-for-cell equality incl. NA positions (the FBM decode returns
  # doubles where the dense path returns integers; values are identical)
  expect_identical(is.na(mr), is.na(md))
  expect_equal(unname(mr), unname(md), ignore_attr = TRUE)
  # frequency: dense and FBM also agree cell for cell
  gl2 <- gl.compliance.check(testset.gl[1:40, 1:80], verbose = 0)
  fb2 <- suppressWarnings(gl.gen2fbm(gl2, verbose = 0))
  set.seed(6)
  ffb <- suppressWarnings(gl.impute(fb2, method = "frequency", verbose = 0))
  set.seed(6)
  fd <- gl.impute(gl2, method = "frequency", verbose = 0)
  mff <- as.matrix(ffb); mfd <- as.matrix(fd)
  expect_identical(is.na(mff), is.na(mfd))
  expect_equal(unname(mff), unname(mfd), ignore_attr = TRUE)
})
