# Characterization tests for gl.select.shapes
# Baseline snapshotted before review (review-gl.select.shapes).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("default returns 0:25 silently", {
  pdf(NULL); on.exit(dev.off())
  o <- capture.output(v <- withVisible(gl.select.shapes(verbose = 0)))
  expect_length(o, 0)
  expect_true(v$visible)   # the vector is the product
  expect_equal(as.numeric(v$value), 0:25)
})

test_that("valid select returned as given; bad ranges caught", {
  pdf(NULL); on.exit(dev.off())
  s <- gl.select.shapes(select = c(1, 1, 5, 8), verbose = 0)
  expect_equal(as.numeric(s), c(1, 1, 5, 8))
  expect_error(gl.select.shapes(select = c(-1, -2), verbose = 0))
  expect_error(gl.select.shapes(select = c(30), verbose = 0))
})

test_that("partially-negative select is caught", {
  # [approved diff F1] baseline: c(-1, 5) slipped the broken
  # min(select < 0 | max(select > 25)) validation and reached pch.
  pdf(NULL); on.exit(dev.off())
  expect_error(gl.select.shapes(select = c(-1, 5), verbose = 0),
               "range 0-25")
})

test_that("x= argument behaviour", {
  # [approved diff F2/F4] baseline: nPop was computed and discarded -
  # wrong-length select passed silently, NULL select returned all 26
  # shapes, and the datatype banner printed at verbose 0.
  pdf(NULL); on.exit(dev.off())
  gl3 <- testset.gl
  pop(gl3) <- factor(rep(c("A", "B", "C"), length.out = nInd(gl3)))
  expect_error(gl.select.shapes(x = gl3, select = c(1, 2), verbose = 0),
               "number of populations")
  o1 <- capture.output(
    s1 <- gl.select.shapes(x = gl3, select = c(1, 2, 5), verbose = 0))
  expect_length(o1, 0)                       # banner silenced at verbose 0
  expect_equal(as.numeric(s1), c(1, 2, 5))
  o2 <- capture.output(s2 <- gl.select.shapes(x = gl3, verbose = 0))
  expect_length(o2, 0)
  expect_equal(as.numeric(s2), 0:2)          # one shape per population
  # more populations than the 26 available shapes -> fatal with guidance
  expect_error(gl.select.shapes(x = testset.gl, verbose = 0),
               "26 shapes")
})

test_that("plot.display suppresses the palette chart", {
  # [approved diff F3] baseline: plot() ran unconditionally.
  grDevices::graphics.off()
  s <- gl.select.shapes(select = c(1, 5), plot.display = FALSE, verbose = 2)
  expect_equal(grDevices::dev.cur(), c("null device" = 1L))
  expect_equal(as.numeric(s), c(1, 5))
})

test_that("input untouched", {
  pdf(NULL); on.exit(dev.off())
  xcopy <- testset.gl
  gl3 <- xcopy
  pop(gl3) <- factor(rep(c("A", "B", "C"), length.out = nInd(gl3)))
  invisible(capture.output(gl.select.shapes(x = gl3, verbose = 0)))
  expect_identical(xcopy, testset.gl)
})
