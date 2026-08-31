# Characterization tests for gl.report.reproducibility
# Baseline snapshotted before review (review-gl.report.reproducibility).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.report.reproducibility returns the input unaltered, invisibly, no history", {
  h <- length(testset.gl@other$history)
  out <- capture.output(
    res <- gl.report.reproducibility(testset.gl, plot.display = FALSE,
                                     verbose = 0)
  )
  expect_identical(res, testset.gl)
  expect_equal(length(res@other$history), h)
  v <- withVisible(gl.report.reproducibility(testset.gl, plot.display = FALSE,
                                             verbose = 0))
  expect_false(v$visible)
})

test_that("gl.report.reproducibility quantile table matches independent recomputation", {
  rep <- testset.gl@other$loc.metrics$RepAvg
  q <- quantile(rep, probs = seq(0, 1, 1 / 20), type = 1, na.rm = TRUE)
  retained <- sapply(q, function(y) sum(rep >= y, na.rm = TRUE))
  out <- capture.output(
    gl.report.reproducibility(testset.gl, plot.display = FALSE, verbose = 1)
  )
  row0 <- out[grepl("\\s0%", out)]
  expect_true(grepl(paste0("\\s", unname(retained["0%"]), "\\s"), row0))
})

test_that("gl.report.reproducibility works on SilicoDArT (Reproducibility metric)", {
  out <- capture.output(
    res <- gl.report.reproducibility(testset.gs, plot.display = FALSE,
                                     verbose = 0)
  )
  expect_identical(res, testset.gs)
})

test_that("gl.report.reproducibility errors when the metric is absent", {
  x <- testset.gl
  x@other$loc.metrics$RepAvg <- NULL
  expect_error(
    gl.report.reproducibility(x, plot.display = FALSE, verbose = 0),
    "RepAvg"
  )
})

test_that("gl.report.reproducibility output at verbose = 0", {
  # [approved diff F2/F3] Pre-fix, 33 lines printed at verbose = 0. Now
  # fully silent.
  out <- capture.output(
    res <- gl.report.reproducibility(testset.gl, plot.display = FALSE,
                                     verbose = 0)
  )
  expect_length(out, 0)
})

test_that("gl.report.reproducibility with plot.file but plot.display = FALSE", {
  # [approved diff F1] Pre-fix this crashed ("object 'p3' not found").
  # Post-fix the RDS save works without displaying.
  tmpdir <- tempdir()
  out <- capture.output(
    res <- gl.report.reproducibility(testset.gl, plot.display = FALSE,
                                     plot.dir = tmpdir,
                                     plot.file = "chartest_rep", verbose = 0)
  )
  expect_identical(res, testset.gl)
  saved <- file.path(tmpdir, "chartest_rep.RDS")
  expect_true(file.exists(saved))
  unlink(saved)
})

test_that("gl.report.reproducibility retained counts are NA-correct", {
  # [approved diff F4] retained counts now exclude NA metric values.
  x <- testset.gl
  x@other$loc.metrics$RepAvg[1:10] <- NA
  out <- capture.output(
    res <- gl.report.reproducibility(x, plot.display = FALSE, verbose = 1)
  )
  row0 <- out[grepl("\\s0%", out)]
  expect_true(grepl("\\s245\\s", row0))
  expect_false(grepl("\\s255\\s", row0))
})
