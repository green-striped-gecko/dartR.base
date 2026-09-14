# Characterization tests for utils.recalc.maf
# Baseline snapshotted before review (recalc-battery review, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.

test_that("maf matches an independent recomputation and verbose 0 is silent", {
  sub <- gl.drop.pop(testset.gl, pop.list = popNames(testset.gl)[1:10],
                     verbose = 0)
  gm <- as.matrix(sub)
  o <- capture.output(r <- utils.recalc.maf(sub, verbose = 0))
  expect_equal(length(o), 0)
  expect_equal(unname(r@other$loc.metrics$maf), unname(pmin(colMeans(gm, na.rm = TRUE) / 2, 1 - colMeans(gm, na.rm = TRUE) / 2)))
  expect_true(r@other$loc.metrics.flags$maf)
})

test_that("an object without a flags slot is handled", {
  sub <- gl.drop.pop(testset.gl, pop.list = popNames(testset.gl)[1:10],
                     verbose = 0)
  sub@other$loc.metrics.flags <- NULL
  # [approved diff B2] baseline: the monomorphs check used `flag == FALSE`,
  # which is `if (logical(0))` when the slot is absent.
  expect_error(capture.output(utils.recalc.maf(sub, verbose = 0)),
    NA)  # [approved diff B2]
})

test_that("SilicoDArT input is rejected", {
  # [approved diff B3] baseline: presence/absence data ran through the
  # diploid arithmetic silently (maf is meaningless for presence/absence data; the doc already says SNP only).
  expect_error(capture.output(utils.recalc.maf(testset.gs, verbose = 0)),
    "SNP")  # [approved diff B3]
})

