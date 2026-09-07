# Characterization tests for gl.reassign.pop
# Baseline snapshotted before review (review-gl.reassign.pop).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("population reassigned from an individual metric", {
  o <- capture.output(g <- gl.reassign.pop(testset.gl, as.pop = "sex",
                                           verbose = 0))
  expect_length(o, 0)
  expect_setequal(popNames(g), c("Female", "Male", "Unknown"))
  expect_equal(nInd(g), nInd(testset.gl))
  expect_equal(length(g@other$history),
               length(testset.gl@other$history) + 1)
  # input untouched
  expect_false("Female" %in% popNames(testset.gl))
})

test_that("nonexistent metric", {
  # [approved diff F1] baseline: no validation - pop(x) was assigned
  # NULL, silently destroying every population assignment. Now fatal.
  expect_error(capture.output(gl.reassign.pop(testset.gl, as.pop = "bogus",
                                              verbose = 0)),
               "not found")
})

test_that("SilicoDArT data accepted", {
  invisible(capture.output(g <- gl.reassign.pop(testset.gs, as.pop = "sex",
                                                verbose = 0)))
  expect_gt(nPop(g), 0)
})
