# Characterization tests for gl.fst.pop
# Baseline snapshotted before review (review-gl.fst.pop), on dartR.data
# 1.2.5 and StAMPP 1.6.3, then updated for the findings approved on
# 2026-09-08. Assertions marked [approved Fn] moved with an approved
# finding recorded in function-review/reports/dartR.base/gl.fst.pop.md;
# every other assertion pins behaviour that must not move.

# ---------------------------------------------------------------- helpers

two_pop <- function() {
  x <- possums.gl[pop(possums.gl) %in% c("A", "B"), ]
  pop(x) <- droplevels(pop(x))
  x
}

# an exchangeable split of one population: true Fst is zero
exchangeable <- function() {
  a <- possums.gl[pop(possums.gl) == "A", ]
  set.seed(7)
  idx <- sample(nInd(a))
  pp <- rep("A2", nInd(a))
  pp[idx[seq_len(nInd(a) %/% 2)]] <- "A1"
  pop(a) <- factor(pp)
  a
}

# raw genlight built from a supplied genotype matrix, no compliance check
raw_gl <- function(m, pv) {
  g <- new("genlight", gen = m,
           ind.names = paste0("i", seq_len(nrow(m))),
           loc.names = paste0("L", seq_len(ncol(m))),
           pop = factor(pv), ploidy = rep(2, nrow(m)))
  g <- methods::as(g, "dartR")
  g@other$loc.metrics <- data.frame(dummy = seq_len(ncol(m)))
  g
}

# ------------------------------------------------------- return structure

test_that("nboots = 1 returns a lower-triangular matrix, not a dist", {
  # [approved F3] @return said "class dist"; the object is a base matrix,
  # and @return now describes it. The class itself is unchanged.
  out <- gl.fst.pop(platypus.gl, nboots = 1, verbose = 0)
  expect_true(is.matrix(out))
  expect_false(inherits(out, "dist"))
  expect_identical(class(out), c("matrix", "array"))
  expect_equal(dim(out), c(3L, 3L))
  expect_true(all(is.na(out[upper.tri(out)])))
  expect_true(all(is.na(diag(out))))
  expect_length(as.numeric(out), 9L)
})

test_that("nboots = 0 is refused", {
  # [approved F6] nboots = 0 previously returned the matrix silently, while
  # @return promised the matrix "if nboots = 1". nboots is now validated as
  # a whole number of 1 or more.
  expect_error(gl.fst.pop(two_pop(), nboots = 0, verbose = 0),
               "nboots must be a single whole number")
})

test_that("nboots > 1 returns the three-element inference list", {
  set.seed(1)
  out <- gl.fst.pop(platypus.gl, nboots = 20, verbose = 0)
  expect_type(out, "list")
  expect_named(out, c("Fsts", "Pvalues", "Bootstraps"))
  expect_true(is.matrix(out$Fsts))
  expect_true(is.matrix(out$Pvalues))
  expect_s3_class(out$Bootstraps, "data.frame")
  expect_equal(dim(out$Bootstraps), c(3L, 26L))
  expect_equal(tail(colnames(out$Bootstraps), 4),
               c("Lower bound CI limit", "Upper bound CI limit",
                 "p-value", "Fst"))
})

# ------------------------------------------------------- pinned estimates

test_that("pairwise theta on platypus.gl is pinned", {
  out <- gl.fst.pop(platypus.gl, nboots = 1, verbose = 0)
  # dimnames follow order of appearance in the data, not factor levels
  expect_identical(rownames(out),
                   c("TENTERFIELD", "SEVERN_BELOW", "SEVERN_ABOVE"))
  expect_identical(rownames(out), colnames(out))
  expect_equal(out[lower.tri(out)],
               c(0.0823972148, 0.0745767976, 0.0600481604),
               tolerance = 1e-8)
})

test_that("pairwise theta on possums.gl is pinned", {
  out <- gl.fst.pop(possums.gl, nboots = 1, verbose = 0)
  expect_identical(rownames(out), LETTERS[1:10])
  expect_equal(sum(out[lower.tri(out)]), 13.4991617506, tolerance = 1e-8)
  expect_equal(out[2, 1], 0.2757253000, tolerance = 1e-8)
  expect_equal(out[10, 1], 0.2916881324, tolerance = 1e-8)
})

test_that("pairwise theta on testset.gl is pinned", {
  out <- gl.fst.pop(testset.gl, nboots = 1, verbose = 0)
  expect_equal(dim(out), c(30L, 30L))
  expect_equal(sum(out[lower.tri(out)], na.rm = TRUE), 176.3633320934,
               tolerance = 1e-8)
})

# ------------------------------------------- independent verification (S1)

test_that("the delegation matches StAMPP::stamppFst called directly", {
  gx <- platypus.gl
  class(gx) <- "genlight"
  direct <- StAMPP::stamppFst(gx, nboots = 1, percent = 95, nclusters = 1)
  out <- gl.fst.pop(platypus.gl, nboots = 1, verbose = 0)
  expect_identical(dimnames(out), dimnames(direct))
  expect_equal(max(abs(out - direct), na.rm = TRUE), 0)
})

test_that("theta matches hierfstat::pairwise.WCfst to machine precision", {
  skip_if_not_installed("hierfstat")
  pl <- gl.filter.monomorphs(
    gl.filter.callrate(platypus.gl, threshold = 1, verbose = 0), verbose = 0)
  out <- gl.fst.pop(pl, nboots = 1, verbose = 0)
  gi <- suppressWarnings(gl2gi(pl, verbose = 0))
  wc <- suppressWarnings(
    hierfstat::pairwise.WCfst(hierfstat::genind2hierfstat(gi),
                              diploid = TRUE))
  pn <- rownames(out)
  for (i in 2:length(pn)) {
    for (j in 1:(i - 1)) {
      expect_equal(out[i, j], wc[pn[i], pn[j]], tolerance = 1e-12,
                   ignore_attr = TRUE)
    }
  }
})

# ------------------------------------------------------- label integrity

test_that("dimnames track content when pop levels are not alphabetical", {
  x <- possums.gl
  pop(x) <- factor(as.character(pop(x)),
                   levels = rev(sort(unique(as.character(pop(x))))))
  out <- gl.fst.pop(x, nboots = 1, verbose = 0)
  pn <- rownames(out)
  for (i in 2:length(pn)) {
    for (j in 1:(i - 1)) {
      sub <- gl.keep.pop(x, pop.list = c(pn[i], pn[j]), verbose = 0)
      v <- gl.fst.pop(sub, nboots = 1, verbose = 0)
      expect_equal(out[i, j], v[lower.tri(v)][1], tolerance = 1e-8)
    }
  }
})

# --------------------------------------------------------- SilicoDArT (F2)

test_that("SilicoDArT is refused", {
  # [approved F2] utils.check.datatype now carries accept = "SNP".
  # Presence/absence data was scored as if every individual were a
  # homozygote and per-population n was halved, and the number returned
  # (0.3032 to 0.6394 on these four populations) was not an Fst.
  gs <- testset.gs[pop(testset.gs) %in% popNames(testset.gs)[1:4], ]
  pop(gs) <- droplevels(pop(gs))
  expect_error(gl.fst.pop(gs, nboots = 1, verbose = 0),
               "found SilicoDArT expecting SNP")
  expect_error(gl.fst.pop(testset.gs, nboots = 1, verbose = 0),
               "found SilicoDArT expecting SNP")
})

# --------------------------------------------------- bootstrap and p-values

test_that("point estimates are unaffected by the bootstrap", {
  set.seed(3)
  a <- gl.fst.pop(possums.gl, nboots = 1, verbose = 0)
  set.seed(4)
  b <- gl.fst.pop(possums.gl, nboots = 10, verbose = 0)$Fsts
  expect_equal(a, b, tolerance = 1e-12)
})

test_that("bootstrap replicates ARE reproducible under set.seed", {
  # [approved F1] StAMPP ran the locus bootstrap inside foreach %dopar% on
  # a PSOCK cluster of its own making; the worker RNG streams were never
  # seeded from the master session, so set.seed() had no effect on the
  # confidence limits or the p-values. The bootstrap is now drawn in the
  # calling session.
  x <- possums.gl[pop(possums.gl) %in% c("A", "B", "C"), ]
  pop(x) <- droplevels(pop(x))
  set.seed(123)
  r1 <- gl.fst.pop(x, nboots = 50, verbose = 0)
  set.seed(123)
  r2 <- gl.fst.pop(x, nboots = 50, verbose = 0)
  set.seed(124)
  r3 <- gl.fst.pop(x, nboots = 50, verbose = 0)
  expect_identical(r1$Bootstraps, r2$Bootstraps)
  expect_identical(r1$Pvalues, r2$Pvalues)
  expect_false(identical(r1$Bootstraps[, 3:52], r3$Bootstraps[, 3:52]))
  expect_false(identical(r1$Bootstraps[["Lower bound CI limit"]],
                         r3$Bootstraps[["Lower bound CI limit"]]))
  # reproducibility comes from the session RNG, not from a new argument
  expect_false("seed" %in% names(formals(gl.fst.pop)))
})

test_that("replicate resampling draws from the full locus set", {
  # [approved F1] the unit of resampling is the locus and the draw runs
  # over 1:nLoc. The sibling defect in gl.report.fstat (PR #384) could only
  # ever reach the first nInd loci of a pair.
  x <- two_pop()
  seen <- integer(0)
  fspy <- gl.fst.pop
  environment(fspy) <- new.env(parent = environment(gl.fst.pop))
  assign("sample", function(x, size, replace) {
    d <- base::sample(x, size, replace = replace)
    seen <<- c(seen, d)
    d
  }, envir = environment(fspy))
  set.seed(77)
  invisible(fspy(x, nboots = 50, verbose = 0))
  expect_equal(length(seen), 50L * nLoc(x))
  expect_equal(max(seen), nLoc(x))
  expect_gt(length(unique(seen)), nInd(x))
  expect_equal(length(unique(seen)), nLoc(x))
})

test_that("a replicate drawn as 1:nLoc reproduces the point estimate", {
  # [approved F1] the replicate estimator and the estimator behind the
  # reported value must be the same one. Intercepting sample() lexically
  # makes every replicate the identity draw.
  fid <- gl.fst.pop
  environment(fid) <- new.env(parent = environment(gl.fst.pop))
  assign("sample", function(x, size, replace) x, envir = environment(fid))
  for (fx in list(two_pop(), exchangeable(),
                  gl.keep.pop(possums.gl, pop.list = c("A", "B", "C"),
                              verbose = 0))) {
    r <- fid(fx, nboots = 3, verbose = 0)
    expect_equal(max(abs(as.matrix(r$Bootstraps[, 3:5]) -
                           r$Fsts[lower.tri(r$Fsts)])), 0)
  }
})

test_that("an exchangeable split gives theta near zero and no signal", {
  # [approved F1] known-answer check on two halves of one population.
  set.seed(11)
  r <- gl.fst.pop(exchangeable(), nboots = 200, verbose = 0)
  expect_lt(abs(r$Fsts[2, 1]), 0.02)
  expect_gt(r$Pvalues[2, 1], 0.05)
  # a genuinely differentiated pair keeps its signal
  set.seed(12)
  d <- gl.fst.pop(two_pop(), nboots = 200, verbose = 0)
  expect_gt(d$Fsts[2, 1], 0.2)
  expect_lt(d$Pvalues[2, 1], 0.05)
  expect_lt(d$Bootstraps[1, "Lower bound CI limit"], d$Bootstraps[1, "Fst"])
  expect_gt(d$Bootstraps[1, "Upper bound CI limit"], d$Bootstraps[1, "Fst"])
})

test_that("the p-value is the one-tailed bootstrap fraction at or below zero", {
  # [DEFECT F5] documented as "probability of Fst values to be different
  # from zero"; implemented as mean(bootstrap replicates <= 0).
  set.seed(11)
  r <- gl.fst.pop(exchangeable(), nboots = 200, verbose = 0)
  bs <- as.numeric(r$Bootstraps[1, 3:202])
  expect_equal(r$Pvalues[2, 1], mean(bs <= 0), ignore_attr = TRUE)
  # exchangeable populations: theta is slightly negative and p is 1, not ~1
  # under a two-sided reading but 1 because every replicate is <= 0
  expect_lt(r$Fsts[2, 1], 0)
  expect_equal(r$Pvalues[2, 1], 1, ignore_attr = TRUE)
})

test_that("the lower CI limit is the sample minimum for nboots <= 40", {
  # The percentile index is ceiling(0.025 * nboots), which is 1 for every
  # nboots up to 40; at nboots = 2 the lower and upper limits collide. The
  # indexing is unchanged; [approved F6] adds the warning that says so.
  b <- gl.fst.pop(two_pop(), nboots = 2, verbose = 0)$Bootstraps
  expect_equal(b[1, "Lower bound CI limit"], b[1, "Upper bound CI limit"])
  b20 <- gl.fst.pop(two_pop(), nboots = 20, verbose = 0)$Bootstraps
  expect_equal(b20[1, "Lower bound CI limit"],
               min(as.numeric(b20[1, 3:22])))
  # [approved F6] the degenerate regime is announced at verbose >= 1
  m20 <- capture.output(invisible(
    gl.fst.pop(two_pop(), nboots = 20, verbose = 1)))
  expect_true(any(grepl("smallest of the 20", m20)))
  m40 <- capture.output(invisible(
    gl.fst.pop(two_pop(), nboots = 40, verbose = 1)))
  expect_true(any(grepl("smallest of the 40", m40)))
  m41 <- capture.output(invisible(
    gl.fst.pop(two_pop(), nboots = 41, verbose = 1)))
  expect_false(any(grepl("smallest of the", m41)))
})

# ----------------------------------------------------- NA and degenerate

test_that("single-individual populations return NaN and are named", {
  # [approved F7] testset.gl ships two populations of n = 1; their pair is
  # NaN. It was silent at every verbosity and is now named at verbose >= 1.
  out <- gl.fst.pop(testset.gl, nboots = 1, verbose = 0)
  cell <- c(out["EmmacNormSalt", "EmmacNormLeic"],
            out["EmmacNormLeic", "EmmacNormSalt"])
  expect_true(any(is.nan(cell)))
  expect_equal(sum(is.nan(out[lower.tri(out)])), 1L)
  msgs <- capture.output(invisible(
    gl.fst.pop(testset.gl, nboots = 1, verbose = 1)))
  expect_true(any(grepl("non-finite", msgs)))
  expect_true(any(grepl("EmmacNormLeic vs EmmacNormSalt", msgs)))
  # still silent at verbose 0
  expect_length(capture.output(invisible(
    gl.fst.pop(testset.gl, nboots = 1, verbose = 0))), 0)
})

test_that("all-monomorphic data returns NaN for every pair, and says so", {
  g <- raw_gl(matrix(0L, 30, 20), rep(c("A", "B", "C"), each = 10))
  out <- gl.fst.pop(g, nboots = 1, verbose = 0)
  expect_true(all(is.nan(out[lower.tri(out)])))
  # [approved F7]
  msgs <- capture.output(invisible(gl.fst.pop(g, nboots = 1, verbose = 1)))
  expect_true(any(grepl("3 of 3 population pairs returned a non-finite",
                        msgs)))
})

test_that("entirely-NA loci and all-NA pop x locus cells stay finite", {
  set.seed(21)
  pv <- rep(c("A", "B", "C"), each = 10)
  m1 <- matrix(sample(0:2, 600, TRUE), 30, 20)
  m1[, 1:5] <- NA
  out1 <- gl.fst.pop(raw_gl(m1, pv), nboots = 1, verbose = 0)
  expect_true(all(is.finite(out1[lower.tri(out1)])))
  m2 <- matrix(sample(0:2, 600, TRUE), 30, 20)
  m2[pv == "A", 1:5] <- NA
  out2 <- gl.fst.pop(raw_gl(m2, pv), nboots = 1, verbose = 0)
  expect_true(all(is.finite(out2[lower.tri(out2)])))
})

test_that("a single locus is accepted", {
  set.seed(22)
  g <- raw_gl(matrix(sample(0:2, 30, TRUE), 30, 1),
              rep(c("A", "B", "C"), each = 10))
  out <- gl.fst.pop(g, nboots = 1, verbose = 0)
  expect_equal(dim(out), c(3L, 3L))
  expect_true(all(is.finite(out[lower.tri(out)])))
})

test_that("two populations with identical genotypes give a negative theta", {
  x <- possums.gl[c(1:5, 1:5), ]
  pop(x) <- factor(rep(c("P", "Q"), each = 5))
  out <- gl.fst.pop(x, nboots = 1, verbose = 0)
  expect_lt(out[2, 1], 0)
})

# --------------------------------------------------------- error paths (F6)

test_that("a single population fails with a dartR message", {
  # [approved F6] was "subscript out of bounds", raised from inside StAMPP.
  one <- gl.keep.pop(possums.gl, pop.list = "A", verbose = 0)
  expect_error(gl.fst.pop(one, nboots = 1, verbose = 0),
               "requires at least two populations")
})

test_that("out-of-range parameters fail with dartR messages", {
  # [approved F6] the FS5 validation block replaces R-internals errors:
  # "subscript out of bounds", "non-numeric argument to binary operator"
  # and "invalid 'length' argument" respectively.
  two <- two_pop()
  expect_error(gl.fst.pop(two, nboots = 10, percent = 150, verbose = 0),
               "percent must be a single number")
  expect_error(gl.fst.pop(two, nboots = 10, percent = "abc", verbose = 0),
               "percent must be a single number")
  expect_error(gl.fst.pop(two, nboots = -5, verbose = 0),
               "nboots must be a single whole number")
  expect_error(gl.fst.pop(two, nboots = 2.5, verbose = 0),
               "nboots must be a single whole number")
  expect_error(gl.fst.pop(two, nclusters = 0, verbose = 0),
               "nclusters must be a single whole number")
})

test_that("non-genlight input is rejected by utils.check.datatype", {
  expect_error(gl.fst.pop(data.frame(a = 1), verbose = 0),
               "inappropriate object passed to function")
})

# --------------------------------------------------------------- contract

test_that("the input object comes back untouched", {
  xcopy <- possums.gl
  invisible(gl.fst.pop(possums.gl, nboots = 1, verbose = 0))
  expect_identical(xcopy, possums.gl)
  expect_identical(class(possums.gl)[1], "dartR")
})

test_that("no history is appended (nothing modified is returned)", {
  nh <- length(possums.gl@other$history)
  invisible(gl.fst.pop(possums.gl, nboots = 1, verbose = 0))
  expect_equal(length(possums.gl@other$history), nh)
})

test_that("verbose = 0 is fully silent", {
  expect_length(
    capture.output(invisible(
      gl.fst.pop(possums.gl, nboots = 1, verbose = 0))), 0)
})

test_that("verbose 3 adds the promised results summary (VRB1)", {
  # [approved F11] levels 2 and 3 were identical at three lines, and only
  # the "Build = Jody" banner distinguished 5. [approved F12] the build tag
  # is gone, so 3 and 5 now agree.
  o1 <- capture.output(invisible(
    gl.fst.pop(possums.gl, nboots = 1, verbose = 1)))
  o2 <- capture.output(invisible(
    gl.fst.pop(possums.gl, nboots = 1, verbose = 2)))
  o3 <- capture.output(invisible(
    gl.fst.pop(possums.gl, nboots = 1, verbose = 3)))
  o5 <- capture.output(invisible(
    gl.fst.pop(possums.gl, nboots = 1, verbose = 5)))
  expect_length(o1, 2L)
  expect_length(o2, 3L)
  expect_gt(length(o3), length(o2))
  expect_true(any(grepl("Populations compared: 10", o3)))
  expect_true(any(grepl("Pairwise comparisons: 45", o3)))
  expect_true(any(grepl("Weir & Cockerham theta", o3)))
  expect_true(any(grepl("Non-finite pairs: 0", o3)))
  expect_identical(o3, o5)
  expect_false(any(grepl("Build = Jody", o5)))
})

test_that("nclusters is honoured and gives the same answer", {
  a <- gl.fst.pop(two_pop(), nboots = 1, nclusters = 1, verbose = 0)
  b <- gl.fst.pop(two_pop(), nboots = 1, nclusters = 2, verbose = 0)
  expect_equal(a, b, tolerance = 1e-12)
})

test_that("the FBM path matches the dense path", {
  skip_if_not(exists("gl.gen2fbm"), "gl.gen2fbm not available")
  x <- gl.keep.pop(possums.gl, pop.list = c("A", "B", "C"), verbose = 0)
  fb <- try(gl.gen2fbm(x, verbose = 0), silent = TRUE)
  skip_if(inherits(fb, "try-error"), "gl.gen2fbm failed on this fixture")
  expect_equal(gl.fst.pop(fb, nboots = 1, verbose = 0),
               gl.fst.pop(x, nboots = 1, verbose = 0), tolerance = 1e-12)
})

# ------------------------------------- downstream consumption (gl.check.panel)

test_that("the gl.check.panel consumption pattern is stable", {
  # dartR.popgen::gl.check.panel:60,62 sorts individuals by population,
  # calls gl.fst.pop(nboots = 1) on the full and panel objects, flattens
  # both with as.numeric() and pairs them with complete.cases(). The pin
  # guards the two properties that pattern depends on: matching dimnames
  # and a matching NA pattern.
  set.seed(31)
  xo <- possums.gl[order(pop(possums.gl)), ]
  xp <- possums.gl[, sample(nLoc(possums.gl), 60)]
  xp <- xp[order(pop(xp)), ]
  fo <- gl.fst.pop(xo, verbose = 0, nboots = 1)
  fp <- gl.fst.pop(xp, verbose = 0, nboots = 1)
  expect_identical(dimnames(fo), dimnames(fp))
  expect_identical(is.na(as.numeric(fo)), is.na(as.numeric(fp)))
  d <- data.frame(o = as.numeric(fo), p = as.numeric(fp))
  d <- d[complete.cases(d), ]
  expect_equal(nrow(d), nPop(xo) * (nPop(xo) - 1) / 2)
  # a real dist return would yield the same 45 paired values (F3 fix is safe)
  expect_setequal(round(d$o, 10), round(as.numeric(as.dist(fo)), 10))
})
