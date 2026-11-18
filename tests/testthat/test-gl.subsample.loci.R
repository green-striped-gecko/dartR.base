# Characterization baseline (function-review, 2026-09-24): outputs of the
# reviewed state, bugs included. Detects change; does not assert
# correctness.
test_that("the deprecated wrapper reproduces gl.subsample.loc", {
  x <- testset.gl
  expect_warning(capture.output(
    v <- gl.subsample.loci(x, n = 50, method = "pic", verbose = 0)), "deprecated")
  top <- order(-x@other$loc.metrics$AvgPIC)[1:50]
  expect_setequal(locNames(v), locNames(x)[top])
  expect_equal(as.character(v@other$loc.metrics$AlleleID),
               as.character(x@other$loc.metrics$AlleleID[
                 match(locNames(v), locNames(x))]))
  # [approved diff, change 5] "PIC" in capitals was sampled at random
  capture.output(v2 <- suppressWarnings(
    gl.subsample.loci(x, n = 50, method = "PIC", verbose = 0)))
  expect_setequal(locNames(v2), locNames(v))
  # the history records the wrapper's call, and replays
  h <- v@other$history[[length(v@other$history)]]
  expect_equal(as.character(h[[1]]), "gl.subsample.loci")
  expect_length(v@other$history, length(x@other$history) + 1)
  capture.output(v3 <- suppressWarnings(eval(h)))
  expect_equal(locNames(v3), locNames(v))
  # out-of-range n still stops, as before
  expect_error(capture.output(suppressWarnings(
    gl.subsample.loci(x, n = 0, verbose = 0))), "subsample size")
  expect_error(capture.output(suppressWarnings(
    gl.subsample.loci(x, n = 1000, verbose = 0))), "subsample size")
  # random sampling without replacement
  capture.output(v4 <- suppressWarnings(gl.subsample.loci(x, n = 30, verbose = 0)))
  expect_equal(nLoc(v4), 30)
  expect_equal(anyDuplicated(locNames(v4)), 0)
})

test_that("a plain adegenet genlight no longer crashes", {
  g <- new("genlight", as.matrix(testset.gl)[1:20, 1:30], ploidy = 2)
  capture.output(v <- suppressWarnings(gl.subsample.loci(g, n = 10, verbose = 0)))
  expect_equal(nLoc(v), 10)
})
