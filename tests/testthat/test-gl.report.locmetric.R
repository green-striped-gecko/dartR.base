# Characterization tests for gl.report.locmetric
# Baseline snapshotted before review (review-gl.report.locmetric).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.report.locmetric returns the input unaltered, invisibly, no history", {
  h <- length(testset.gl@other$history)
  out <- capture.output(
    res <- gl.report.locmetric(testset.gl, metric = "SnpPosition",
                               plot.display = FALSE, verbose = 0)
  )
  expect_identical(res, testset.gl)
  expect_equal(length(res@other$history), h)
  v <- withVisible(gl.report.locmetric(testset.gl, metric = "SnpPosition",
                                       plot.display = FALSE, verbose = 0))
  expect_false(v$visible)
})

test_that("gl.report.locmetric quantile table matches independent recomputation", {
  v <- testset.gl@other$loc.metrics$SnpPosition
  q <- quantile(v, probs = seq(0, 1, 1 / 20), type = 1, na.rm = TRUE)
  retained <- sapply(q, function(y) sum(v >= y, na.rm = TRUE))
  out <- capture.output(
    gl.report.locmetric(testset.gl, metric = "SnpPosition",
                        plot.display = FALSE, verbose = 1)
  )
  row0 <- out[grepl("\\s0%", out)]
  expect_true(grepl(paste0("\\s", unname(retained["0%"]), "\\s"), row0))
})

test_that("gl.report.locmetric works on SilicoDArT and a custom metric", {
  out <- capture.output(
    res <- gl.report.locmetric(testset.gs, metric = "AvgReadDepth",
                               plot.display = FALSE, verbose = 0)
  )
  expect_identical(res, testset.gs)
})

test_that("gl.report.locmetric error paths", {
  expect_error(
    gl.report.locmetric(testset.gl, metric = "nope", verbose = 0),
    "not found"
  )
  expect_error(
    gl.report.locmetric(testset.gl, metric = "AlleleID", verbose = 0),
    "not numeric"
  )
})

test_that("gl.report.locmetric plot.file works without display", {
  tmpdir <- tempdir()
  out <- capture.output(
    res <- gl.report.locmetric(testset.gl, metric = "SnpPosition",
                               plot.display = FALSE, plot.dir = tmpdir,
                               plot.file = "chartest_lm", verbose = 0)
  )
  saved <- file.path(tmpdir, "chartest_lm.RDS")
  expect_true(file.exists(saved))
  unlink(saved)
})

test_that("gl.report.locmetric output at verbose = 0", {
  # [approved diff F1] Pre-fix, 32 lines printed at verbose = 0. Now fully
  # silent.
  out <- capture.output(
    res <- gl.report.locmetric(testset.gl, metric = "SnpPosition",
                               plot.display = FALSE, verbose = 0)
  )
  expect_length(out, 0)
})

test_that("gl.report.locmetric stats lines formatting", {
  # [approved diff F2] summary() now runs on the vector; stats lines carry
  # a single clean label and numeric value.
  out <- capture.output(
    gl.report.locmetric(testset.gl, metric = "SnpPosition",
                        plot.display = FALSE, verbose = 1)
  )
  expect_false(any(grepl("Min\\.", out)))
  expect_true(any(grepl("Minimum\\s+:\\s+5", out)))
  expect_true(any(grepl("3rd quartile", out)))
})

test_that("gl.report.locmetric retained counts are NA-correct", {
  # [approved diff F3] retained counts now exclude NA metric values.
  x <- testset.gl
  x@other$loc.metrics$SnpPosition[1:10] <- NA
  out <- capture.output(
    res <- gl.report.locmetric(x, metric = "SnpPosition",
                               plot.display = FALSE, verbose = 1)
  )
  row0 <- out[grepl("\\s0%", out)]
  expect_true(grepl("\\s245\\s", row0))
  expect_false(grepl("\\s255\\s", row0))
})
