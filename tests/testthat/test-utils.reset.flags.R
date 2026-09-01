# Characterization tests for utils.reset.flags
# Baseline snapshotted before review (recalc-battery review, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.

test_that("flags are reset and verbose 0 is silent", {
  sub <- gl.drop.pop(testset.gl, pop.list = popNames(testset.gl)[1:10],
                     verbose = 0)
  o <- capture.output(r <- utils.reset.flags(sub, set = FALSE, verbose = 0))
  expect_equal(length(o), 0)
  expect_false(r@other$loc.metrics.flags$AvgPIC)
  expect_false(r@other$loc.metrics.flags$maf)
  o <- capture.output(r2 <- utils.reset.flags(sub, set = TRUE, verbose = 0))
  expect_true(r2@other$loc.metrics.flags$CallRate)
})

test_that("SilicoDArT branch resets its own flags and disables SNP flags", {
  o <- capture.output(r <- utils.reset.flags(testset.gs, set = TRUE,
                                             verbose = 0))
  expect_true(r@other$loc.metrics.flags$OneRatio)
  expect_false(r@other$loc.metrics.flags$PICRef)
})

test_that("the value argument sets the stored verbosity", {
  sub <- gl.drop.pop(testset.gl, pop.list = popNames(testset.gl)[1:10],
                     verbose = 0)
  sub@other$verbose <- NULL
  # [approved diff R1] baseline: value was validated then ignored; the
  # verbosity slot was hardcoded to 2.
  o <- capture.output(r <- utils.reset.flags(sub, value = 4, verbose = 0))
  expect_equal(r@other$verbose, 4)  # [approved diff R1]
})

test_that("an out-of-range value warns only at verbose >= 1", {
  sub <- gl.drop.pop(testset.gl, pop.list = popNames(testset.gl)[1:10],
                     verbose = 0)
  # [approved diff R2] baseline: the coercion warning printed at verbose 0.
  o <- capture.output(r <- utils.reset.flags(sub, value = 9, verbose = 0))
  expect_equal(length(o), 0)  # [approved diff R2]
})

test_that("no bogus loc.metrics$monomorphs column is invented", {
  sub <- gl.drop.pop(testset.gl, pop.list = popNames(testset.gl)[1:10],
                     verbose = 0)
  sub@other$loc.metrics$monomorphs <- NULL
  # [approved diff R3] baseline: the SNP branch created a loc.metrics
  # column for what is a flag, not a locus metric.
  o <- capture.output(r <- utils.reset.flags(sub, verbose = 0))
  expect_false("monomorphs" %in% names(r@other$loc.metrics))  # [approved diff R3]
})
