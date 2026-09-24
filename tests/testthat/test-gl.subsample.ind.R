# Characterization baseline (function-review, 2026-09-24): outputs of the
# reviewed state, bugs included. Detects change; does not assert
# correctness.
test_that("baseline: sample sizes returned", {
  x <- testset.gl
  y <- platypus.gl  # populations of 23, 17 and 41
  size <- function(...) {
    capture.output(v <- suppressWarnings(gl.subsample.ind(..., verbose = 0)))
    v
  }
  # [approved diff, change 3] the default n = NULL crashed; it is now the
  # smallest population size (1 in testset.gl) per population
  expect_equal(nInd(size(x)), nPop(x))
  expect_equal(nInd(size(x, n = 30, by.pop = FALSE, replace = TRUE)), 30)
  # [approved diff, change 2] was 255 (nLoc) and an error
  expect_equal(nInd(size(x, n = 300, by.pop = FALSE, replace = TRUE)), 300)
  expect_equal(nInd(size(x, n = 600, by.pop = FALSE, replace = TRUE)), 600)
  expect_equal(nInd(size(x, n = 252, by.pop = FALSE, replace = FALSE)), 250)
  expect_equal(as.vector(table(pop(size(y, n = 10, by.pop = TRUE)))),
               c(10, 10, 10))
  # [approved diff, change 1] was 93/87/70
  expect_equal(as.vector(table(pop(size(y, n = 70, by.pop = TRUE)))),
               c(70, 70, 70))
  expect_equal(as.vector(table(pop(size(y, n = 100, by.pop = TRUE)))),
               c(100, 100, 100))
  v <- size(y, n = 10, by.pop = TRUE, replace = FALSE)
  expect_equal(as.vector(table(pop(v))), c(10, 10, 10))
  expect_equal(anyDuplicated(indNames(v)), 0)
  expect_equal(nrow(v@other$ind.metrics), nInd(v))
  expect_length(v@other$history, length(y@other$history) + 1)
})

test_that("invalid n stops with a message naming it", {
  expect_error(capture.output(gl.subsample.ind(testset.gl, n = 0, verbose = 0)),
               "n, the number of individuals")
  expect_error(capture.output(gl.subsample.ind(testset.gl, n = "a", verbose = 0)),
               "n, the number of individuals")
})
