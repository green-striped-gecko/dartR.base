# Characterization tests for gl.report.taglength
# Baseline snapshotted before review (review-gl.report.taglength).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.report.taglength returns the input unaltered, invisibly, no history", {
  h <- length(testset.gl@other$history)
  out <- capture.output(
    res <- gl.report.taglength(testset.gl, plot.display = FALSE, verbose = 0)
  )
  expect_identical(res, testset.gl)
  expect_equal(length(res@other$history), h)
  v <- withVisible(gl.report.taglength(testset.gl, plot.display = FALSE,
                                       verbose = 0))
  expect_false(v$visible)
})

test_that("gl.report.taglength quantile table matches independent recomputation", {
  ncht <- nchar(as.character(testset.gl@other$loc.metrics$TrimmedSequence))
  q <- quantile(ncht, probs = seq(0, 1, 1 / 20), type = 1, na.rm = TRUE)
  retained <- sapply(q, function(y) sum(ncht >= y, na.rm = TRUE))
  out <- capture.output(
    gl.report.taglength(testset.gl, plot.display = FALSE, verbose = 1)
  )
  row0 <- out[grepl("\\s0%", out)]
  expect_true(grepl(paste0("\\s", unname(retained["0%"]), "\\s"), row0))
  row100 <- out[grepl("\\s100%", out)]
  expect_true(grepl(paste0("\\s", unname(retained["100%"]), "\\s"), row100))
})

test_that("gl.report.taglength works on SilicoDArT with TrimmedSequence", {
  out <- capture.output(
    res <- gl.report.taglength(testset.gs, plot.display = FALSE, verbose = 0)
  )
  expect_identical(res, testset.gs)
})

test_that("gl.report.taglength errors when TrimmedSequence is absent", {
  x <- testset.gl
  x@other$loc.metrics$TrimmedSequence <- NULL
  expect_error(
    gl.report.taglength(x, plot.display = FALSE, verbose = 0),
    "TrimmedSequence"
  )
})

test_that("gl.report.taglength output at verbose = 0", {
  # [approved diff F2/F3] Pre-fix, 33 lines printed at verbose = 0. Now
  # fully silent (assigned inside capture.output; return is invisible).
  out <- capture.output(
    res <- gl.report.taglength(testset.gl, plot.display = FALSE, verbose = 0)
  )
  expect_length(out, 0)
})

test_that("gl.report.taglength with plot.file but plot.display = FALSE", {
  # [approved diff F1] Pre-fix this crashed ("object 'p3' not found").
  # Post-fix the RDS save works without displaying.
  tmpdir <- tempdir()
  out <- capture.output(
    res <- gl.report.taglength(testset.gl, plot.display = FALSE,
                               plot.dir = tmpdir,
                               plot.file = "chartest_taglen", verbose = 0)
  )
  expect_identical(res, testset.gl)
  saved <- file.path(tmpdir, "chartest_taglen.RDS")
  expect_true(file.exists(saved))
  unlink(saved)
})

test_that("gl.report.taglength retained counts are NA-correct", {
  # [approved diff F4] retained counts now exclude NA tag lengths.
  x <- testset.gl
  x@other$loc.metrics$TrimmedSequence[1:10] <- NA
  out <- capture.output(
    res <- gl.report.taglength(x, plot.display = FALSE, verbose = 1)
  )
  row0 <- out[grepl("\\s0%", out)]
  expect_true(grepl("\\s245\\s", row0))
  expect_false(grepl("\\s255\\s", row0))
})
