# Characterization tests for gl.sim.genotypes (dartR.base)
#
# Captured during the function-review campaign at upstream/dev = ddaed27, then
# updated when findings F1-F14 were approved on 2026-09-09. Blocks that changed
# because an approved finding changed behaviour are marked "[approved Fn]" at
# the expectation that moved.
#
# Report: function-review/reports/dartR.base/gl.sim.genotypes.md

library(dartR.data)

# ---- fixtures ---------------------------------------------------------------

# possums.gl: 300 individuals, 200 loci, 10 populations, no all-NA locus.
# n.ind is now per population, so a call on this object returns 10 * n.ind
# individuals.
sg_src <- possums.gl

# A single population, for the cases where pooling and structure must not be
# confounded.
sg_pop1 <- gl.keep.pop(possums.gl, pop.list = popNames(possums.gl)[1],
                       verbose = 0)

# A wide single population (1000 loci), so that the model probes can use
# n.ind = 1000 without the n.ind <= nLoc cap (F1) binding.
sg_wide <- gl.keep.pop(platypus.gl, pop.list = popNames(platypus.gl)[1],
                       verbose = 0)

# Build a small genlight straight from a dosage matrix.
sg_fixture <- function(mat) {
  rownames(mat) <- paste0("I", seq_len(nrow(mat)))
  colnames(mat) <- paste0("L", seq_len(ncol(mat)))
  g <- new("genlight",
           gen = mat,
           ind.names = rownames(mat),
           loc.names = colnames(mat),
           ploidy = rep(2, nrow(mat)))
  gl.compliance.check(g, verbose = 0)
}

# ---- returned object --------------------------------------------------------

test_that("returns a dartR/genlight of the requested shape", {
  set.seed(1)
  out <- gl.sim.genotypes(sg_src, n.ind = 50, verbose = 0)

  expect_s4_class(out, "genlight")
  # [approved F3] n.ind is per population: 50 individuals in each of 10.
  expect_equal(nInd(out), 50 * nPop(sg_src))
  expect_equal(nLoc(out), nLoc(sg_src))
  expect_identical(locNames(out), locNames(sg_src))
  expect_identical(indNames(out), paste0("Ind_", 1:(50 * nPop(sg_src))))
  expect_equal(length(unique(indNames(out))), 50 * nPop(sg_src))

  # DAT1: ploidy is 2 for every individual and the vector is per-individual.
  expect_equal(length(ploidy(out)), nInd(out))
  expect_true(all(ploidy(out) == 2))

  # SNP dosages only. possums.gl has no all-NA locus in any population, so no
  # missing data is simulated from it.
  expect_setequal(sort(unique(as.vector(as.matrix(out)))), c(0, 1, 2))
  expect_equal(sum(is.na(as.matrix(out))), 0)
})

test_that("metadata is constructed and tracks the genotypes 1:1", {
  set.seed(1)
  out <- gl.sim.genotypes(sg_src, n.ind = 20, verbose = 0)

  expect_equal(nrow(out@other$loc.metrics), nLoc(out))
  expect_equal(nrow(out@other$ind.metrics), nInd(out))

  # [approved F9] ind.metrics carries id and pop, so a downstream rebuild of
  # pop from ind.metrics has something to read.
  expect_identical(names(out@other$ind.metrics), c("id", "pop"))
  expect_identical(out@other$ind.metrics$id, indNames(out))
  expect_identical(out@other$ind.metrics$pop, as.character(pop(out)))

  # [approved F9] no column or flag literally named "array(NA, nLoc(x))" /
  # "array(NA, 1)": loc.metrics and loc.metrics.flags are supplied shaped, so
  # gl.compliance.check has nothing to invent.
  expect_false("array(NA, nLoc(x))" %in% names(out@other$loc.metrics))
  expect_false("array(NA, 1)" %in% names(out@other$loc.metrics.flags))
  expect_true(all(c("CallRate", "maf") %in% names(out@other$loc.metrics)))

  # Derived metrics describe the SIMULATED data, not the source: callrate is 1
  # because no missing data is generated from a source with no all-NA locus.
  expect_true(all(out@other$loc.metrics$CallRate == 1))
})

test_that("history records the call that created the object", {
  # [approved F5] FS8: the returned object carries gl.sim.genotypes(...) as its
  # only history entry. It previously carried only the internal
  # gl.recalc.metrics and gl.compliance.check calls.
  set.seed(1)
  out <- gl.sim.genotypes(sg_src, n.ind = 20, verbose = 0)
  h <- vapply(out@other$history, function(e) as.character(e[[1]]), character(1))

  expect_equal(length(out@other$history), 1L)
  expect_identical(h, "gl.sim.genotypes")
  expect_true(grepl("n.ind = 20", paste(deparse(out@other$history[[1]]),
                                        collapse = " ")))
})

test_that("source population structure is carried through", {
  # [approved F3] the source has 10 populations; the result has the same 10,
  # under the same names, with n.ind individuals in each.
  set.seed(1)
  out <- gl.sim.genotypes(sg_src, n.ind = 20, verbose = 0)
  expect_equal(nPop(sg_src), 10L)
  expect_identical(popNames(out), popNames(sg_src))
  expect_true(all(table(pop(out)) == 20))

  set.seed(2)
  out1 <- gl.sim.genotypes(sg_pop1, n.ind = 20, verbose = 0)
  expect_identical(popNames(out1), popNames(sg_pop1))
})

# ---- the simulation model ---------------------------------------------------

test_that("genotypes are drawn under Hardy-Weinberg from the locus frequencies", {
  set.seed(99)
  out <- gl.sim.genotypes(sg_wide, n.ind = 1000, verbose = 0)
  M <- as.matrix(out)
  f_in <- gl.allele.freq(sg_wide, by = "loc", verbose = 0)$frequency
  f_sim <- colMeans(M, na.rm = TRUE) / 2

  ok <- !is.na(f_in) & !is.na(f_sim)
  expect_gt(sum(ok), 900)

  # Allele frequencies are recovered, with the same polarity as the source
  # (the simulated dosage counts the same allele gl.allele.freq reports).
  expect_gt(cor(f_in[ok], f_sim[ok]), 0.99)
  expect_lt(max(abs(f_in[ok] - f_sim[ok])), 0.05)

  # Per-locus HWE, 1 df, against the simulated frequency.
  chi <- vapply(seq_len(ncol(M)), function(j) {
    p <- f_sim[j]
    col <- M[, j]
    if (is.na(p) || p <= 0 || p >= 1 || all(is.na(col))) return(NA_real_)
    col <- col[!is.na(col)]
    o <- tabulate(col + 1, 3)
    e <- length(col) * c((1 - p)^2, 2 * p * (1 - p), p^2)
    sum((o - e)^2 / e)
  }, numeric(1))
  chi <- chi[!is.na(chi)]
  expect_gt(length(chi), 300)
  expect_lt(mean(chi), 1.6)          # expectation 1
  expect_lt(mean(chi > 3.841), 0.12) # expectation 0.05
})

test_that("JOINT STRUCTURE: loci with equal allele frequency are independent", {
  # Recycling-class probe. A short random vector recycled across the matrix
  # leaves every per-locus marginal correct but makes congruent columns
  # identical. Correlations between equal-frequency locus pairs must sit near 0.
  set.seed(2024)
  M <- as.matrix(gl.sim.genotypes(sg_wide, n.ind = 1000, verbose = 0))
  f <- round(gl.allele.freq(sg_wide, by = "loc", verbose = 0)$frequency, 4)

  shared <- as.numeric(names(table(f))[table(f) >= 2])
  cors <- unlist(lapply(shared, function(v) {
    idx <- which(f == v)
    if (length(idx) < 2) return(NULL)
    unlist(lapply(seq_len(length(idx) - 1), function(a)
      vapply((a + 1):length(idx), function(b)
        suppressWarnings(cor(M[, idx[a]], M[, idx[b]], use = "complete.obs")),
        numeric(1))))
  }))
  cors <- cors[!is.na(cors)]

  expect_gt(length(cors), 50)
  expect_lt(mean(abs(cors)), 0.10)
  expect_lt(max(abs(cors)), 0.30)
})

test_that("JOINT STRUCTURE: unrelated individuals match only at the chance rate", {
  set.seed(2025)
  M <- as.matrix(gl.sim.genotypes(sg_wide, n.ind = 1000, verbose = 0))
  M <- M[, colSums(is.na(M)) == 0, drop = FALSE]

  # Chance rate = mean over loci of sum_g p_g^2, from the simulated marginals.
  expected <- mean(apply(M, 2, function(col) {
    p <- tabulate(col + 1, 3) / length(col)
    sum(p^2)
  }))
  observed <- replicate(500, {
    i <- sample(nrow(M), 1)
    j <- sample(setdiff(seq_len(nrow(M)), i), 1)
    mean(M[i, ] == M[j, ])
  })

  expect_lt(abs(mean(observed) - expected), 0.02)
})

test_that("JOINT STRUCTURE: two independent draws are consumed per genotype", {
  # The RNG stream must advance by exactly 2 draws per locus per population,
  # each of length n.ind: 2 draws per simulated cell and no reuse.
  # [approved F3] the draw count is now per population, so it scales with nPop.
  n.ind <- 40
  n.loc <- nLoc(sg_src)
  n.pop <- nPop(sg_src)

  set.seed(7)
  invisible(gl.sim.genotypes(sg_src, n.ind = n.ind, verbose = 0))
  after_fn <- .Random.seed

  set.seed(7)
  for (k in seq_len(2 * n.loc * n.pop)) {
    invisible(sample(c(0, 1), size = n.ind, replace = TRUE, prob = c(0.4, 0.6)))
  }
  after_manual <- .Random.seed

  expect_identical(after_fn, after_manual)

  # And the result equals an explicit two-haplotype reimplementation, drawn
  # from each population's own allele frequencies.
  m <- gl.allele.freq(sg_src, by = "popxloc", verbose = 0)
  m$popn <- as.character(m$popn)
  set.seed(7)
  v <- matrix(NA_real_, nrow = n.loc, ncol = n.ind * n.pop)
  for (k in seq_len(n.pop)) {
    mp <- m[m$popn == popNames(sg_src)[k], ]
    f <- mp[order(mp$loc_order), "frequency"]
    cols <- ((k - 1) * n.ind + 1):(k * n.ind)
    for (i in seq_len(n.loc)) {
      if (is.na(f[i])) next
      v1 <- sample(c(0, 1), size = n.ind, replace = TRUE,
                   prob = c(1 - f[i], f[i]))
      v2 <- sample(c(0, 1), size = n.ind, replace = TRUE,
                   prob = c(1 - f[i], f[i]))
      v[i, cols] <- v1 + v2
    }
  }
  set.seed(7)
  out <- gl.sim.genotypes(sg_src, n.ind = n.ind, verbose = 0)
  expect_equal(unname(t(v)),
               unname(matrix(as.matrix(out), nrow = n.ind * n.pop)))
})

test_that("allele frequencies are estimated within each population", {
  # [approved F3] the description says the genotypes are simulated from the
  # allele frequencies of a population. Each population is now simulated from
  # its own frequencies, so simulated Ho tracks each SOURCE population's Ho
  # instead of collapsing to the pooled He (0.471).
  Msrc <- as.matrix(sg_src)
  set.seed(11)
  sim <- gl.sim.genotypes(sg_src, n.ind = 200, verbose = 0)
  Msim <- as.matrix(sim)

  Ho_src <- vapply(popNames(sg_src), function(p)
    mean(colMeans(Msrc[pop(sg_src) == p, , drop = FALSE] == 1, na.rm = TRUE),
         na.rm = TRUE), numeric(1))
  Ho_sim <- vapply(popNames(sg_src), function(p)
    mean(colMeans(Msim[pop(sim) == p, , drop = FALSE] == 1, na.rm = TRUE),
         na.rm = TRUE), numeric(1))

  # Per-population tracking: every population lands within 0.05 of its source.
  expect_true(all(abs(Ho_sim - Ho_src) < 0.05))
  expect_gt(cor(Ho_src, Ho_sim), 0.95)

  # The overall mean no longer inflates to the pooled He.
  f_pool <- colMeans(Msrc, na.rm = TRUE) / 2
  He_pooled <- mean(2 * f_pool * (1 - f_pool))
  expect_gt(He_pooled, 0.45)                       # 0.471
  expect_lt(abs(mean(Ho_sim) - mean(Ho_src)), 0.05)
  expect_lt(mean(Ho_sim), 0.40)

  # A genuinely single-population source is reproduced faithfully.
  Ho_p1 <- mean(colMeans(as.matrix(sg_pop1) == 1, na.rm = TRUE))
  set.seed(5)
  Ho_p1sim <- mean(colMeans(as.matrix(
    gl.sim.genotypes(sg_pop1, n.ind = 200, verbose = 0)) == 1))
  expect_lt(abs(Ho_p1sim - Ho_p1), 0.05)
})

test_that("monomorphic loci stay monomorphic", {
  # 20 individuals x 60 loci: 20 fixed at 0, 20 fixed at 2, 20 at frequency 0.5.
  # Wider than n.ind so the F1 cap does not bind.
  m <- cbind(matrix(0, 20, 20), matrix(2, 20, 20), matrix(1, 20, 20))
  g <- sg_fixture(m)
  expect_equal(unique(gl.allele.freq(g, by = "loc", verbose = 0)$frequency),
               c(0, 1, 0.5))

  set.seed(2)
  out <- gl.sim.genotypes(g, n.ind = 50, verbose = 0)
  M <- as.matrix(out)
  expect_identical(unique(as.vector(M[, 1:20])), 0L)
  expect_identical(unique(as.vector(M[, 21:40])), 2L)
  expect_setequal(sort(unique(as.vector(M[, 41:60]))), c(0L, 1L, 2L))
})

# ---- the n.ind parameter ----------------------------------------------------

test_that("the n.ind > n.loc cap is applied", {
  # [approved F1] the function prints "Setting n.ind to <n.loc>" and now
  # simulates n.loc individuals per population; previously it printed the
  # message and simulated the full requested n.ind.
  n.loc <- nLoc(sg_pop1)
  set.seed(1)
  msg <- capture.output(
    out <- gl.sim.genotypes(sg_pop1, n.ind = n.loc + 100, verbose = 1))

  expect_true(any(grepl("Setting n.ind to", msg)))
  expect_true(any(grepl(paste("Setting n.ind to", n.loc), msg)))
  expect_equal(nInd(out), n.loc)

  # On a structured source the cap applies to each population.
  set.seed(1)
  outm <- gl.sim.genotypes(sg_src, n.ind = nLoc(sg_src) + 100, verbose = 0)
  expect_equal(nInd(outm), nLoc(sg_src) * nPop(sg_src))
})

test_that("n.ind edge values", {
  set.seed(1)
  one <- gl.sim.genotypes(sg_pop1, n.ind = 1, verbose = 0)
  expect_equal(nInd(one), 1)
  expect_equal(length(ploidy(one)), 1)

  # [approved F7] invalid values are named and rejected by the function, not
  # by base R ("subscript out of bounds", "negative length vectors").
  expect_error(gl.sim.genotypes(sg_pop1, n.ind = 0, verbose = 0),
               "n.ind must be 1 or greater")
  expect_error(gl.sim.genotypes(sg_pop1, n.ind = -5, verbose = 0),
               "n.ind must be 1 or greater")
  expect_error(gl.sim.genotypes(sg_pop1, n.ind = "a", verbose = 0),
               "n.ind must be a single finite numeric value")
  expect_error(gl.sim.genotypes(sg_pop1, n.ind = NA, verbose = 0),
               "n.ind must be a single finite numeric value")
  expect_error(gl.sim.genotypes(sg_pop1, n.ind = c(10, 20), verbose = 0),
               "n.ind must be a single finite numeric value")
  expect_error(gl.sim.genotypes(sg_pop1, n.ind = NULL, verbose = 0),
               "n.ind must be a single finite numeric value")
  expect_error(gl.sim.genotypes(sg_pop1, n.ind = Inf, verbose = 0),
               "n.ind must be a single finite numeric value")

  # [approved F7] a non-integer n.ind is rounded, with a gated warning, not
  # silently truncated by array().
  set.seed(1)
  msg <- capture.output(
    frac <- gl.sim.genotypes(sg_pop1, n.ind = 2.7, verbose = 1))
  expect_equal(nInd(frac), 3)
  expect_true(any(grepl("Rounding", msg)))
})

# ---- input handling ---------------------------------------------------------

test_that("an all-NA locus is simulated as an all-NA locus", {
  # [approved F2] gl.allele.freq returns NA for a locus with no calls in a
  # population. That locus is now left missing for the individuals of that
  # population instead of aborting inside sample(). testset.gl, testset2.gl and
  # platypus.gl all carry such loci, so the function previously failed on four
  # of the five packaged datasets.
  for (d in list(testset.gl, testset2.gl, platypus.gl)) {
    set.seed(1)
    out <- gl.sim.genotypes(d, n.ind = 5, verbose = 0)
    expect_s4_class(out, "genlight")

    # nLoc is preserved and the all-NA loci of the source are the all-NA loci
    # of the result.
    expect_equal(nLoc(out), nLoc(d))
    src_allna <- which(colSums(!is.na(as.matrix(d))) == 0)
    sim_allna <- which(colSums(!is.na(as.matrix(out))) == 0)
    expect_gt(length(src_allna), 0)
    expect_identical(as.integer(sim_allna), as.integer(src_allna))
  }

  # Missingness is reproduced population by population: a locus is missing in
  # the simulated individuals of population p exactly when it has no calls in
  # population p of the source. Nothing else is simulated as missing.
  set.seed(1)
  out <- gl.sim.genotypes(testset.gl, n.ind = 5, verbose = 0)
  Msrc <- as.matrix(testset.gl)
  Msim <- as.matrix(out)
  src_gap <- vapply(popNames(testset.gl), function(p)
    colSums(!is.na(Msrc[pop(testset.gl) == p, , drop = FALSE])) == 0,
    logical(nLoc(testset.gl)))
  sim_gap <- vapply(popNames(out), function(p)
    colSums(!is.na(Msim[pop(out) == p, , drop = FALSE])) == 0,
    logical(nLoc(out)))
  expect_identical(sim_gap, src_gap)
  expect_equal(sum(is.na(Msim)), sum(src_gap) * 5)

  # A source with no within-population gaps simulates no missing data at all.
  set.seed(1)
  out <- gl.sim.genotypes(sg_wide, n.ind = 20, verbose = 0)
  expect_equal(sum(is.na(as.matrix(out))),
               sum(colSums(!is.na(as.matrix(sg_wide))) == 0) * 20)
})

test_that("SilicoDArT input is refused", {
  # [approved F4] DAT7: presence/absence data has no diploid dosage. The
  # function previously read the P/A rate as an allele frequency and returned a
  # ploidy-2 0/1/2 object that utils.check.datatype reported as SNP.
  expect_equal(utils.check.datatype(testset.gs, verbose = 0), "SilicoDArT")
  expect_error(gl.sim.genotypes(testset.gs, n.ind = 20, verbose = 0),
               "found SilicoDArT expecting SNP")
})

test_that("a plain genlight is accepted", {
  # [approved F10] DAT5: the input is passed through gl.compliance.check, so a
  # genlight not built by dartR no longer dies inside gl.allele.freq with
  # "argument is of length zero".
  raw <- new("genlight",
             gen = as.matrix(sg_src)[1:10, 1:10],
             ploidy = rep(2, 10))
  set.seed(1)
  out <- gl.sim.genotypes(raw, n.ind = 5, verbose = 0)
  expect_s4_class(out, "genlight")
  expect_equal(nInd(out), 5)
  expect_equal(nLoc(out), 10)
})

test_that("the input object is not modified", {
  before <- sg_src
  set.seed(1)
  invisible(gl.sim.genotypes(sg_src, n.ind = 10, verbose = 0))
  expect_identical(before, sg_src)
})

# ---- reproducibility and verbosity ------------------------------------------

test_that("the same seed gives the same object", {
  set.seed(123); a <- gl.sim.genotypes(sg_src, n.ind = 30, verbose = 0)
  set.seed(123); b <- gl.sim.genotypes(sg_src, n.ind = 30, verbose = 0)
  set.seed(124); c <- gl.sim.genotypes(sg_src, n.ind = 30, verbose = 0)

  expect_identical(as.matrix(a), as.matrix(b))
  expect_identical(a@other$loc.metrics, b@other$loc.metrics)
  expect_false(identical(as.matrix(a), as.matrix(c)))
})

test_that("verbose = 0 is silent for an ordinary call", {
  set.seed(1)
  expect_silent(invisible(gl.sim.genotypes(sg_src, n.ind = 10, verbose = 0)))
})

test_that("verbose = 0 is silent when the n.ind warning fires", {
  # [approved F6] VRB3/VRB5: the warning was an ungated cat(), so verbose = 0
  # printed 2 lines.
  set.seed(1)
  out <- capture.output(
    invisible(gl.sim.genotypes(sg_pop1, n.ind = nLoc(sg_pop1) + 1,
                               verbose = 0)))
  expect_equal(length(out), 0)
})

test_that("verbose levels 1 to 5 print the expected banners", {
  set.seed(1)
  v1 <- capture.output(invisible(gl.sim.genotypes(sg_pop1, n.ind = 10, verbose = 1)))
  v2 <- capture.output(invisible(gl.sim.genotypes(sg_pop1, n.ind = 10, verbose = 2)))
  v3 <- capture.output(invisible(gl.sim.genotypes(sg_pop1, n.ind = 10, verbose = 3)))
  v5 <- capture.output(invisible(gl.sim.genotypes(sg_pop1, n.ind = 10, verbose = 5)))

  expect_equal(length(v1), 2)
  expect_true(any(grepl("Starting gl.sim.genotypes", v1)))
  expect_true(any(grepl("Completed: gl.sim.genotypes", v1)))
  expect_equal(length(v2), 3)

  # [approved F12] verbose = 3 now prints the results summary VRB1 promises;
  # its output was previously identical to verbose = 2.
  expect_false(identical(v2, v3))
  expect_true(any(grepl("Simulated 10 individuals", v3)))
  expect_true(any(grepl("Mean simulated heterozygosity", v3)))

  # [approved F13] FS3: the outdated build argument is gone, so verbose = 5 no
  # longer reports "Build = v.2023.3" on a package at a different version.
  expect_false(any(grepl("Build = v.2023.3", v5)))
})

# ---- downstream usability ---------------------------------------------------

test_that("the returned object survives compliance and downstream functions", {
  set.seed(1)
  out <- gl.sim.genotypes(sg_src, n.ind = 30, verbose = 0)

  expect_s4_class(gl.compliance.check(out, verbose = 0), "genlight")
  expect_s4_class(gl.filter.callrate(out, threshold = 0.9, verbose = 0),
                  "genlight")
  expect_s4_class(gl.recalc.metrics(out, verbose = 0), "genlight")

  # A simulated object carrying all-NA loci also survives call-rate filtering.
  set.seed(1)
  na_out <- gl.sim.genotypes(testset.gl, n.ind = 10, verbose = 0)
  expect_s4_class(gl.filter.callrate(na_out, threshold = 0.9, verbose = 0),
                  "genlight")
})
