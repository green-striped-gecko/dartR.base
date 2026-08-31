# Characterization tests for gl.report.overshoot
# Baseline snapshotted before review (review-gl.report.overshoot).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.report.overshoot count matches independent recomputation", {
  trimmed <- as.character(testset.gl@other$loc.metrics$TrimmedSequence)
  snpos <- testset.gl@other$loc.metrics$SnpPosition
  expected <- sum((snpos + 1) > nchar(trimmed))
  expect_equal(expected, 21)  # testset.gl has 21 genuine overshoot loci
  out <- capture.output(res <- gl.report.overshoot(testset.gl, verbose = 1))
  expect_true(any(grepl(paste0("outside the trimmed sequence: ", expected), out)))
})

test_that("gl.report.overshoot detects crafted overshoot loci", {
  xo <- testset.gl
  xo@other$loc.metrics$SnpPosition[1:3] <- 1000
  out <- capture.output(res <- gl.report.overshoot(xo, verbose = 1))
  expect_true(any(grepl("outside the trimmed sequence: 24", out)))
  expect_true(any(grepl(locNames(testset.gl)[1], out, fixed = TRUE)))
})

test_that("gl.report.overshoot returns the input unaltered, invisibly, no history", {
  h <- length(testset.gl@other$history)
  out <- capture.output(res <- gl.report.overshoot(testset.gl, verbose = 0))
  expect_identical(res, testset.gl)
  expect_equal(length(res@other$history), h)
  v <- withVisible(gl.report.overshoot(testset.gl, verbose = 0))
  expect_false(v$visible)
})

test_that("gl.report.overshoot rejects SilicoDArT and missing metrics", {
  expect_error(gl.report.overshoot(testset.gs, verbose = 0))
  x6 <- testset.gl
  x6@other$loc.metrics$TrimmedSequence <- NULL
  expect_error(gl.report.overshoot(x6, verbose = 0), "TrimmedSequence")
  x7 <- testset.gl
  x7@other$loc.metrics$SnpPosition <- NULL
  expect_error(gl.report.overshoot(x7, verbose = 0), "position")
})

test_that("gl.report.overshoot output at verbose = 0", {
  # [approved diff F1] Pre-fix, the results printed unconditionally
  # (2 lines at verbose = 0). Now fully silent.
  out <- capture.output(res <- gl.report.overshoot(testset.gl, verbose = 0))
  expect_length(out, 0)
})

test_that("gl.report.overshoot listing has no stray trailing comma", {
  # [approved diff F2] listing now uses paste(collapse = ", ").
  out <- capture.output(res <- gl.report.overshoot(testset.gl, verbose = 1))
  listing <- out[grepl("100050384", out)]
  expect_false(grepl(",$", listing))
  expect_true(grepl("100114332-43-A/G$", listing))
})
