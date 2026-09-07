# Characterization tests for gl.alf
# Baseline snapshotted before review (allele-frequency family redundancy
# review, dev at ddaed27). Assertions tagged [approved Fn] were flipped
# when the matching approved finding was applied; every other assertion is
# unchanged from the baseline and still passes.

test_that("gl.alf returns a two-column data.frame keyed by locus name", {
  a <- gl.alf(testset.gl)
  expect_s3_class(a, "data.frame")
  expect_equal(dim(a), c(nLoc(testset.gl), 2L))
  expect_equal(colnames(a), c("alf1", "alf2"))
  # [approved F3] @return previously named the columns "ref" and "alt";
  # it now names the columns the function actually produces, alf1/alf2
  expect_identical(rownames(a), locNames(testset.gl))
  expect_true(all(vapply(a, is.numeric, logical(1))))
})

test_that("gl.alf equals hand computation colMeans(as.matrix(x))/2 (SNP)", {
  a <- gl.alf(testset.gl)
  m <- as.matrix(testset.gl)
  h <- unname(colMeans(m, na.rm = TRUE) / 2)
  expect_equal(a$alf2, h)
  expect_equal(a$alf1, 1 - h)
  # alf1 + alf2 is exactly 1 wherever defined
  expect_true(all(abs(a$alf1 + a$alf2 - 1) < 1e-12, na.rm = TRUE))
})

test_that("pinned values on testset.gl anchor the numbers", {
  a <- gl.alf(testset.gl)
  expect_equal(a["100049687-12-C/T", "alf2"], 1)
  expect_equal(a["100049698-16-G/A", "alf2"], 0.03947368, tolerance = 1e-7)
  expect_equal(a["100049728-23-A/G", "alf2"], 0)
  expect_equal(range(a$alf2, na.rm = TRUE), c(0, 1))
  expect_equal(sum(a$alf1 == 1, na.rm = TRUE), 79L)
  expect_equal(sum(a$alf1 == 0, na.rm = TRUE), 62L)
})

test_that("pinned values on possums.gl (the @examples dataset)", {
  a <- gl.alf(possums.gl[, 1:10])
  expect_equal(rownames(a), paste0("X", 1:10))
  expect_equal(a$alf2[1], 0.5366667, tolerance = 1e-6)
  expect_equal(a$alf1[6], 0.3416667, tolerance = 1e-6)
})

test_that("all-NA loci yield NaN in both columns", {
  # testset.gl carries three loci with no scored genotypes
  a <- gl.alf(testset.gl)
  expect_equal(sum(is.nan(a$alf1)), 3L)
  expect_equal(sum(is.nan(a$alf2)), 3L)
  expect_identical(which(is.nan(a$alf1)), which(is.nan(a$alf2)))

  # explicit fixture: locus 3 blanked
  m <- as.matrix(testset.gl[, 1:5])
  m[, 3] <- NA
  fx <- new("genlight", m, ploidy = 2)
  r <- gl.alf(fx)
  expect_true(is.nan(r$alf1[3]))
  expect_true(is.nan(r$alf2[3]))
  expect_false(anyNA(r$alf1[-3]))
})

test_that("NA policy is per locus, na.rm = TRUE, over the whole object", {
  # a single individual leaves NaN wherever that individual is unscored
  si <- testset.gl[1, ]
  r <- gl.alf(si)
  m <- as.matrix(si)
  expect_equal(sum(is.nan(r$alf2)), sum(is.na(m[1, ])))
  expect_equal(r$alf2, unname(m[1, ]) / 2)
})

test_that("single-locus input returns a one-row frame", {
  r <- gl.alf(testset.gl[, 1])
  expect_equal(dim(r), c(1L, 2L))
  expect_identical(rownames(r), locNames(testset.gl)[1])
})

test_that("gl.alf accepts a plain genlight not built by dartR", {
  # no preamble, so no loc.metrics.flags access: non-dartR input works
  m <- matrix(sample(0:2, 200, TRUE), nrow = 10)
  colnames(m) <- paste0("L", 1:20)
  plain <- new("genlight", m, ploidy = 2)
  r <- gl.alf(plain)
  expect_s3_class(r, "data.frame")
  expect_identical(rownames(r), paste0("L", 1:20))
})

test_that("gl.alf is silent and returns visibly", {
  expect_silent(invisible(gl.alf(testset.gl[, 1:5])))
  expect_length(capture.output(invisible(gl.alf(testset.gl[, 1:5]))), 0L)
  expect_true(withVisible(gl.alf(testset.gl[, 1:2]))$visible)
})

test_that("duplicate locus names keep the locus keys", {
  # [approved F1] the baseline pinned rownames == as.character(1:20):
  # data.frame() fell back to 1:n whenever the colMeans names were not
  # unique, and callers that read rownames() as locus names
  # (gl.report.heterozygosity, gl.select.panel) got integers. The row keys
  # are now set from locNames() and disambiguated with make.unique().
  d <- testset.gl[, 1:20]
  locNames(d)[5] <- locNames(d)[4]
  r <- gl.alf(d)
  expect_identical(rownames(r), make.unique(locNames(d)))
  expect_identical(rownames(r)[-5], locNames(d)[-5])
  expect_identical(rownames(r)[5], paste0(locNames(d)[5], ".1"))
  expect_false(any(rownames(r) == as.character(seq_len(20))))
  # values are positional and unchanged by the row keys
  expect_equal(r$alf2, unname(colMeans(as.matrix(d), na.rm = TRUE) / 2))
})

test_that("duplicate locus names no longer corrupt gl.report.heterozygosity", {
  # [approved F1] the baseline produced polyLoc = nLoc for every population
  # and negative monoLoc counts (-1, -2) through
  # gl.report.heterozygosity.r:508, which reads rownames(gl.alf(...))
  d <- testset.gl[, 1:30]
  locNames(d)[5] <- locNames(d)[4]
  h <- gl.report.heterozygosity(d, method = "pop", verbose = 0,
                                plot.display = FALSE)
  expect_true(all(h$monoLoc >= 0))
  expect_true(all(h$polyLoc < nLoc(d)))
  expect_true(all(h$monoLoc + h$polyLoc <= nLoc(d)))
})

test_that("SilicoDArT is rejected by the datatype gate", {
  # [approved F2] the baseline accepted Tag P/A data and divided it by the
  # SNP ploidy, returning half the presence frequency (alf2 capped at 0.5,
  # alf1 never below 0.5). It is now a fatal error.
  expect_equal(unique(ploidy(testset.gs)), 1)
  expect_error(gl.alf(testset.gs), "SilicoDArT")
  expect_error(gl.alf(testset.gs), "expecting SNP")
})

test_that("FBM-backed input gives the same answer as the in-memory object", {
  skip_if_not(exists("gl.gen2fbm"), "gl.gen2fbm not available")
  skip_if_not_installed("bigstatsr")
  f <- tryCatch(gl.gen2fbm(possums.gl, verbose = 0), error = function(e) NULL)
  skip_if(is.null(f), "FBM conversion unavailable")
  r <- gl.alf(f)
  expect_equal(r, gl.alf(possums.gl))
  expect_identical(rownames(r), locNames(possums.gl))
})

# ---------------------------------------------------------------------
# Cross-function equivalence: gl.alf vs gl.allele.freq(simple = TRUE)
# These pin the CURRENT relationship between the two functions. The
# SilicoDArT assertion is expected to flip when gl.allele.freq PR #374
# (F2: SilicoDArT divisor in the by = 'loc' overwrite) is applied.
# ---------------------------------------------------------------------

test_that("SNP: the two agree to 4 dp but are not identical (rounding)", {
  a <- gl.alf(testset.gl)
  s <- gl.allele.freq(testset.gl, simple = TRUE, verbose = 0)
  expect_identical(dim(a), dim(s))
  expect_identical(colnames(a), colnames(s))
  expect_identical(rownames(a), rownames(s))
  expect_false(identical(a$alf2, s$alf2))
  # gl.allele.freq rounds to 4 dp; gl.alf does not
  expect_equal(s$alf2, round(a$alf2, 4))
  expect_lte(max(abs(a$alf2 - s$alf2), na.rm = TRUE), 5e-05)
  expect_gt(sum(a$alf2 != s$alf2, na.rm = TRUE), 0)
  # NaN loci coincide
  expect_identical(which(is.nan(a$alf2)), which(is.nan(s$alf2)))
})

test_that("SNP: gl.allele.freq(by = 'loc') frequency equals gl.alf alf2", {
  a <- gl.alf(testset.gl)
  l <- gl.allele.freq(testset.gl, by = "loc", verbose = 0)
  expect_equal(l$frequency, round(a$alf2, 4))
})

test_that("SilicoDArT: gl.alf refuses, gl.allele.freq still answers", {
  # [approved F2] the baseline pinned the two agreeing at half frequency on
  # ploidy-1 data. gl.alf now refuses Tag P/A input outright, so the two
  # functions no longer share a SilicoDArT domain. gl.allele.freq's own
  # divisor is the subject of its PR #374 and is not asserted here.
  expect_error(gl.alf(testset.gs))
  s <- gl.allele.freq(testset.gs, simple = TRUE, verbose = 0)
  expect_s3_class(s, "data.frame")
})

test_that("gl.allele.freq(simple = TRUE) rejects input gl.alf accepts", {
  # non-dartR genlight: gl.alf works, gl.allele.freq errors on the
  # unguarded loc.metrics.flags access (gl.allele.freq F3)
  m <- matrix(sample(0:2, 200, TRUE), nrow = 10)
  colnames(m) <- paste0("L", 1:20)
  plain <- new("genlight", m, ploidy = 2)
  expect_s3_class(gl.alf(plain), "data.frame")
  expect_error(gl.allele.freq(plain, simple = TRUE, verbose = 0))

  # duplicate locus names: gl.alf degrades, gl.allele.freq errors
  d <- testset.gl[, 1:20]
  locNames(d)[5] <- locNames(d)[4]
  expect_s3_class(gl.alf(d), "data.frame")
  expect_error(gl.allele.freq(d, simple = TRUE, verbose = 0))
})

# ---------------------------------------------------------------------
# Cross-function equivalence: the utils.recalc.freq* decomposition
# ---------------------------------------------------------------------

test_that("FreqHomSnp + FreqHets/2 reproduces gl.alf alf2 exactly", {
  y <- gl.compliance.check(testset.gl, verbose = 0)
  y <- utils.recalc.freqhomsnp(y, verbose = 0)
  y <- utils.recalc.freqhomref(y, verbose = 0)
  y <- utils.recalc.freqhets(y, verbose = 0)
  lm <- y@other$loc.metrics
  expect_true(all(abs(lm$FreqHomRef + lm$FreqHets + lm$FreqHomSnp - 1) < 1e-12,
                  na.rm = TRUE))
  expect_equal(lm$FreqHomSnp + lm$FreqHets / 2, gl.alf(y)$alf2)
})

test_that("utils.recalc.maf is a pure transform of gl.alf alf2", {
  y <- gl.compliance.check(testset.gl, verbose = 0)
  a <- gl.alf(y)$alf2
  y2 <- utils.recalc.maf(y, verbose = 0)
  mf <- y2@other$loc.metrics$maf
  expect_equal(mf, ifelse(a > 0.5, 1 - a, a))
  # gl.alf emits NaN for all-NA loci; ifelse converts them to NA
  expect_equal(sum(is.na(mf)), 3L)
  expect_equal(sum(is.nan(mf)), 0L)
})

test_that("per-population gl.alf matches gl.allele.freq popxloc to 4 dp", {
  sub <- gl.keep.pop(testset.gl, popNames(testset.gl)[1:3], verbose = 0)
  sp <- seppop(sub)
  alf_pop <- sapply(sp, function(p) gl.alf(p)$alf2)
  f <- gl.allele.freq(sub, percent = TRUE, by = "popxloc", verbose = 0)
  fm <- matrix(f$frequency, nrow = nLoc(sub), byrow = TRUE)
  expect_lte(max(abs(as.vector(alf_pop) - as.vector(fm) / 100), na.rm = TRUE),
             1e-04)
  # both drop missing genotypes within a cell; empty cells are NaN/NA
  expect_equal(sum(is.nan(alf_pop)), sum(is.na(fm)))
})

test_that("gl.tree.nj's own frequency step uses a different NA policy", {
  # [pins defect in gl.tree.nj, reported there] mean(e)/2 without
  # na.rm propagates NA far more widely than the gl.alf policy
  sub <- gl.keep.pop(testset.gl, popNames(testset.gl)[1:3], verbose = 0)
  tn <- apply(as.matrix(sub), 2, tapply, pop(sub), function(e) mean(e) / 2)
  sp <- seppop(sub)
  alf_pop <- sapply(sp, function(p) gl.alf(p)$alf2)
  expect_gt(sum(is.na(tn)), sum(is.nan(alf_pop)))
})
