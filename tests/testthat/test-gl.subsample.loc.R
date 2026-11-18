# Characterization baseline (function-review, 2026-09-24): outputs of the
# reviewed state. Detects change; does not assert correctness.
test_that("baseline: loci sampled and metrics kept in step", {
  x <- testset.gl
  set.seed(1)
  capture.output(v <- gl.subsample.loc(x, n = 50, replace = FALSE, verbose = 0))
  expect_equal(nLoc(v), 50)
  expect_equal(as.character(v@other$loc.metrics$AlleleID),
               as.character(x@other$loc.metrics$AlleleID[
                 match(locNames(v), locNames(x))]))
  capture.output(v <- gl.subsample.loc(x, n = 255, replace = TRUE, verbose = 0))
  expect_equal(nLoc(v), 255)
  expect_equal(anyDuplicated(locNames(v)), 0)
  expect_equal(nrow(v@other$loc.metrics), 255)
  # n above nLoc is capped
  capture.output(v <- suppressWarnings(gl.subsample.loc(x, n = 500, verbose = 0)))
  expect_equal(nLoc(v), 255)
  expect_error(capture.output(gl.subsample.loc(x, verbose = 0)),
               "n, the number of loci to subsample, must be supplied")
})

test_that("method = 'pic' picks the most informative loci, recalculated when stale", {
  x <- testset.gl
  capture.output(v <- gl.subsample.loc(x, n = 50, method = "pic", verbose = 0))
  expect_setequal(locNames(v),
                  locNames(x)[order(-x@other$loc.metrics$AvgPIC)[1:50]])
  expect_equal(anyDuplicated(locNames(v)), 0)
  # after dropping populations the stored AvgPIC is out of date
  x2 <- gl.drop.pop(x, popNames(x)[1:20], verbose = 0)
  expect_false(isTRUE(x2@other$loc.metrics.flags$AvgPIC))
  x3 <- gl.recalc.metrics(x2, verbose = 0)
  capture.output(v <- gl.subsample.loc(x2, n = 50, method = "PIC", verbose = 0))
  expect_setequal(locNames(v),
                  locNames(x3)[order(-x3@other$loc.metrics$AvgPIC)[1:50]])
  # SilicoDArT ranks on PIC
  capture.output(v <- gl.subsample.loc(testset.gs, n = 20, method = "pic", verbose = 0))
  expect_equal(nLoc(v), 20)
  expect_error(capture.output(gl.subsample.loc(x, n = 5, method = "bogus", verbose = 0)),
               "method must be")
})

test_that("mono.rm removes monomorphic loci before sampling, history records one call", {
  x <- testset.gl
  capture.output(xm <- gl.filter.monomorphs(x, verbose = 0))
  capture.output(v <- gl.subsample.loc(x, n = nLoc(xm), replace = FALSE,
                                       mono.rm = TRUE, verbose = 0))
  expect_setequal(locNames(v), locNames(xm))
  expect_length(v@other$history, length(x@other$history) + 1)
  expect_equal(as.character(v@other$history[[length(v@other$history)]][[1]]),
               "gl.subsample.loc")
})
