# Characterization baseline (function-review, 2026-09-24): outputs of the
# reviewed state. Detects change; does not assert correctness.
test_that("baseline: significant tests per population on bandicoot.gl", {
  skip_if_not_installed("HardyWeinberg")
  capture.output(o <- gl.hwe.pop(bandicoot.gl, plot.out = FALSE, verbose = 0))
  expect_equal(dim(o$HWE), c(5, 1000))
  expect_equal(unname(rowSums(o$HWE)), c(22, 27, 53, 25, 20))
  expect_null(o$plot)
  expect_null(o$HWformat)
})

test_that("baseline: HWformat returns genotype counts per population", {
  skip_if_not_installed("HardyWeinberg")
  capture.output(o <- gl.hwe.pop(bandicoot.gl, plot.out = FALSE,
                                 HWformat = TRUE, verbose = 0))
  expect_length(o$HWformat, 5)
  expect_named(o$HWformat[[1]], c("AA", "AB", "BB"))
  expect_equal(nrow(o$HWformat[[1]]), 1000)
})

test_that("an object without populations is tested as one population", {
  skip_if_not_installed("HardyWeinberg")
  x <- bandicoot.gl
  pop(x) <- NULL
  capture.output(o <- gl.hwe.pop(x, plot.out = FALSE, verbose = 0))
  expect_equal(dim(o$HWE), c(1, nLoc(x)))
  expect_equal(rownames(o$HWE), "pop1")
})
