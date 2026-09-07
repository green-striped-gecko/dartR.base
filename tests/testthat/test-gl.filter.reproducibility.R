# Characterization tests for gl.filter.reproducibility
# Baseline snapshotted before review (review-gl.filter.reproducibility).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.filter.reproducibility retains loci with repeatability >= threshold", {
  rep <- testset.gl@other$loc.metrics$RepAvg
  expected <- sum(rep >= 0.99)
  res <- gl.filter.reproducibility(testset.gl, threshold = 0.99,
                                   plot.display = FALSE, verbose = 0)
  expect_equal(nLoc(res), expected)
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
  expect_true(all(res@other$loc.metrics$RepAvg >= 0.99))
  # boundary loci retained
  expect_equal(sum(res@other$loc.metrics$RepAvg == 0.99), sum(rep == 0.99))
})

test_that("gl.filter.reproducibility works on SilicoDArT (Reproducibility)", {
  reps <- testset.gs@other$loc.metrics$Reproducibility
  res <- gl.filter.reproducibility(testset.gs, threshold = 0.99,
                                   plot.display = FALSE, verbose = 0)
  expect_equal(nLoc(res), sum(reps >= 0.99))
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
})

test_that("gl.filter.reproducibility leaves the input untouched", {
  x <- testset.gl
  invisible(gl.filter.reproducibility(x, plot.display = FALSE, verbose = 0))
  expect_identical(x, testset.gl)
})

test_that("gl.filter.reproducibility history entries", {
  # [approved diff F1] Pre-fix, the internal gl.drop.loc delegation leaked
  # its own history entry (one call -> two entries). Post-fix exactly one
  # entry is appended, and it is this function's own call.
  h <- length(testset.gl@other$history)
  res <- gl.filter.reproducibility(testset.gl, threshold = 0.99,
                                   plot.display = FALSE, verbose = 0)
  expect_equal(length(res@other$history), h + 1)
  expect_match(deparse(res@other$history[[h + 1]])[1],
               "gl.filter.reproducibility")
})

test_that("gl.filter.reproducibility with the metric absent", {
  # [approved diff F2] Pre-fix this was a SILENT NO-OP (all loci returned
  # unfiltered). Post-fix it stops with the same error as the report
  # sibling.
  xx <- testset.gl
  xx@other$loc.metrics$RepAvg <- NULL
  expect_error(
    gl.filter.reproducibility(xx, plot.display = FALSE, verbose = 0),
    "RepAvg"
  )
})

test_that("gl.filter.reproducibility with NA metric values", {
  # [approved diff F4] Pre-fix, NA-RepAvg loci silently passed the filter.
  # Post-fix they are removed (aligned with gl.filter.rdepth/taglength).
  xna <- testset.gl
  xna@other$loc.metrics$RepAvg[1:10] <- NA
  res <- gl.filter.reproducibility(xna, threshold = 0.99,
                                   plot.display = FALSE, verbose = 0)
  expect_equal(sum(is.na(res@other$loc.metrics$RepAvg)), 0)
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
  rep_clean <- testset.gl@other$loc.metrics$RepAvg
  expect_equal(nLoc(res), sum(rep_clean[-(1:10)] >= 0.99))
})

test_that("gl.filter.reproducibility with plot.file but plot.display = FALSE", {
  # [approved diff F3] Pre-fix this crashed ("object 'p3' not found").
  # Post-fix the RDS save works without displaying.
  tmpdir <- tempdir()
  res <- gl.filter.reproducibility(testset.gl, plot.display = FALSE,
                                   plot.dir = tmpdir,
                                   plot.file = "chartest_frep", verbose = 0)
  saved <- file.path(tmpdir, "chartest_frep.RDS")
  expect_true(file.exists(saved))
  unlink(saved)
})

test_that("gl.filter.reproducibility threshold warning at verbose = 0", {
  # [approved diff F5] warnings now gate at verbose >= 1.
  out <- capture.output(
    res <- gl.filter.reproducibility(testset.gl, threshold = 2,
                                     plot.display = FALSE, verbose = 0)
  )
  expect_length(out, 0)
})

test_that("gl.filter.reproducibility returns invisibly", {
  # [approved diff F6] now returns invisibly, matching the other filters.
  v <- withVisible(gl.filter.reproducibility(testset.gl, plot.display = FALSE,
                                             verbose = 0))
  expect_false(v$visible)
})
