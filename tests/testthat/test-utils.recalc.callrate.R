# Characterization tests for utils.recalc.callrate
# Baseline snapshotted before review (recalc-battery review, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.

test_that("CallRate matches an independent recomputation and verbose 0 is silent", {
  sub <- gl.drop.pop(testset.gl, pop.list = popNames(testset.gl)[1:10],
                     verbose = 0)
  gm <- as.matrix(sub)
  o <- capture.output(r <- utils.recalc.callrate(sub, verbose = 0))
  expect_equal(length(o), 0)
  expect_equal(unname(r@other$loc.metrics$CallRate), unname(signif(1 - colSums(is.na(gm)) / nInd(sub), 6)))
  expect_true(r@other$loc.metrics.flags$CallRate)
})

test_that("an object without a flags slot is handled (isTRUE idiom)", {
  sub <- gl.drop.pop(testset.gl, pop.list = popNames(testset.gl)[1:10],
                     verbose = 0)
  sub@other$loc.metrics.flags <- NULL
  expect_error(capture.output(utils.recalc.callrate(sub, verbose = 0)), NA)
})

test_that("SilicoDArT is supported (datatype-agnostic metric)", {
  o <- capture.output(r <- utils.recalc.callrate(testset.gs, verbose = 0))
  expect_equal(unname(r@other$loc.metrics$CallRate),
    unname(signif(1 - colSums(is.na(as.matrix(testset.gs))) / nInd(testset.gs), 6)))
})
