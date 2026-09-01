# Characterization tests for utils.recalc.freqhomref
# Baseline snapshotted before review (recalc-battery review, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.

test_that("FreqHomRef matches an independent recomputation and verbose 0 is silent", {
  sub <- gl.drop.pop(testset.gl, pop.list = popNames(testset.gl)[1:10],
                     verbose = 0)
  gm <- as.matrix(sub)
  o <- capture.output(r <- utils.recalc.freqhomref(sub, verbose = 0))
  expect_equal(length(o), 0)
  expect_equal(unname(r@other$loc.metrics$FreqHomRef), unname(colMeans(gm == 0, na.rm = TRUE)))
  expect_true(r@other$loc.metrics.flags$FreqHomRef)
})

test_that("an object without a flags slot is handled", {
  sub <- gl.drop.pop(testset.gl, pop.list = popNames(testset.gl)[1:10],
                     verbose = 0)
  sub@other$loc.metrics.flags <- NULL
  # [approved diff B2] baseline: the monomorphs check used `flag == FALSE`,
  # which is `if (logical(0))` when the slot is absent.
  expect_error(capture.output(utils.recalc.freqhomref(sub, verbose = 0)),
    NA)  # [approved diff B2]
})

test_that("SilicoDArT input is rejected", {
  # [approved diff B3] baseline: presence/absence data ran through the
  # diploid arithmetic silently (absences counted as homozygous reference).
  expect_error(capture.output(utils.recalc.freqhomref(testset.gs, verbose = 0)),
    "SNP")  # [approved diff B3]
})

