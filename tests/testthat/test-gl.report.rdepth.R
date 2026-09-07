# Characterization tests for gl.report.rdepth
# Baseline snapshotted before review (review-gl.report.rdepth).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.report.rdepth returns the input unaltered, invisibly, no history", {
  h <- length(testset.gl@other$history)
  out <- capture.output(
    res <- gl.report.rdepth(testset.gl, plot.display = FALSE, verbose = 0)
  )
  expect_identical(res, testset.gl)
  expect_equal(length(res@other$history), h)
  v <- withVisible(gl.report.rdepth(testset.gl, plot.display = FALSE,
                                    verbose = 0))
  expect_false(v$visible)
})

test_that("gl.report.rdepth quantile table matches independent recomputation", {
  rdepth <- testset.gl@other$loc.metrics$rdepth
  q <- quantile(rdepth, probs = seq(0, 1, 1 / 20), type = 1, na.rm = TRUE)
  retained <- sapply(q, function(y) sum(rdepth >= y, na.rm = TRUE))
  # verbose = 1: post-F2 the table prints at verbose >= 1 (silent at 0)
  out <- capture.output(
    gl.report.rdepth(testset.gl, plot.display = FALSE, verbose = 1)
  )
  # the printed table ends with the 0% row: all loci retained
  expect_true(any(grepl("^21\\s+0%", out)))
  row0 <- out[grepl("^21\\s+0%", out)]
  expect_true(grepl(paste0("\\s", unname(retained["0%"]), "\\s"), row0))
  # 100% row retains the loci at the maximum read depth
  row100 <- out[grepl("^1\\s+100%", out)]
  expect_true(grepl(paste0("\\s", unname(retained["100%"]), "\\s"), row100))
})

test_that("gl.report.rdepth works on SilicoDArT (AvgReadDepth)", {
  out <- capture.output(
    res <- gl.report.rdepth(testset.gs, plot.display = FALSE, verbose = 0)
  )
  expect_identical(res, testset.gs)
})

test_that("gl.report.rdepth output at verbose = 0", {
  # [approved diff F2/F3] Pre-fix, the stats block + quantile table printed
  # unconditionally (33 lines at verbose = 0). Now gated at verbose >= 1;
  # verbose = 0 is fully silent (assigned inside capture.output because the
  # return is invisible anyway).
  out <- capture.output(
    res <- gl.report.rdepth(testset.gl, plot.display = FALSE, verbose = 0)
  )
  expect_length(out, 0)
})

test_that("gl.report.rdepth with plot.file but plot.display = FALSE", {
  # [approved diff F1] Pre-fix this crashed ("object 'p3' not found").
  # Post-fix the plots are always built and the RDS save works without
  # displaying.
  tmpdir <- tempdir()
  out <- capture.output(
    res <- gl.report.rdepth(testset.gl, plot.display = FALSE,
                            plot.dir = tmpdir,
                            plot.file = "chartest_rdepth", verbose = 0)
  )
  expect_identical(res, testset.gl)
  saved <- file.path(tmpdir, "chartest_rdepth.RDS")
  expect_true(file.exists(saved))
  unlink(saved)
})

test_that("gl.report.rdepth retained counts are NA-correct", {
  # [approved diff F4] retained counts now exclude NA read depths.
  x <- testset.gl
  x@other$loc.metrics$rdepth[1:10] <- NA
  out <- capture.output(
    res <- gl.report.rdepth(x, plot.display = FALSE, verbose = 1)
  )
  row0 <- out[grepl("\\s0%", out)]
  expect_true(grepl("\\s245\\s", row0))
  expect_false(grepl("\\s255\\s", row0))
})

test_that("gl.report.rdepth errors when the read-depth metric is absent", {
  x <- testset.gl
  x@other$loc.metrics$rdepth <- NULL
  expect_error(
    gl.report.rdepth(x, plot.display = FALSE, verbose = 0),
    "Read depth not included"
  )
})
