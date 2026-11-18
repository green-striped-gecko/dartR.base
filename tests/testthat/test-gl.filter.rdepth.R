# Characterization tests for gl.filter.rdepth
# Baseline snapshotted before review (review-gl.filter.rdepth).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.filter.rdepth retains loci within [lower, upper] inclusive", {
  rdepth <- testset.gl@other$loc.metrics$rdepth
  expected <- sum(rdepth >= 8 & rdepth <= 50)
  res <- gl.filter.rdepth(testset.gl, lower = 8, upper = 50,
                          plot.display = FALSE, verbose = 0)
  expect_equal(nLoc(res), expected)
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
  expect_true(all(res@other$loc.metrics$rdepth >= 8 &
                    res@other$loc.metrics$rdepth <= 50))
  # boundary values are retained (documented: removed only below/above)
  expect_equal(sum(res@other$loc.metrics$rdepth == 8),
               sum(rdepth == 8))
})

test_that("gl.filter.rdepth works on SilicoDArT (AvgReadDepth)", {
  ard <- testset.gs@other$loc.metrics$AvgReadDepth
  res <- gl.filter.rdepth(testset.gs, lower = 8, upper = 50,
                          plot.display = FALSE, verbose = 0)
  expect_equal(nLoc(res), sum(ard >= 8 & ard <= 50))
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
})

test_that("gl.filter.rdepth appends history and leaves the input untouched", {
  x <- testset.gl
  h <- length(x@other$history)
  res <- gl.filter.rdepth(x, plot.display = FALSE, verbose = 0)
  expect_identical(x, testset.gl)
  expect_equal(length(res@other$history), h + 1)
})

test_that("gl.filter.rdepth is silent at verbose = 0", {
  out <- capture.output(
    res <- gl.filter.rdepth(testset.gl, plot.display = FALSE, verbose = 0)
  )
  expect_length(out, 0)
})

test_that("gl.filter.rdepth with NA read depths", {
  # [approved diff F1] Pre-fix, NA rdepth values silently DESYNCED
  # genotypes and metadata (127 loci vs 137 metric rows, 10 garbage all-NA
  # rows). Post-fix NA-depth loci are removed cleanly and sync holds.
  xna <- testset.gl
  xna@other$loc.metrics$rdepth[1:10] <- NA
  res <- gl.filter.rdepth(xna, lower = 8, upper = 50,
                          plot.display = FALSE, verbose = 0)
  expect_equal(nLoc(res), 127)
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
  expect_equal(sum(is.na(res@other$loc.metrics$rdepth)), 0)
})

test_that("gl.filter.rdepth with plot.file but plot.display = FALSE", {
  # [approved diff F2] Pre-fix this crashed ("object 'p3' not found").
  # Post-fix the plots are always built and the RDS save works without
  # displaying.
  tmpdir <- tempdir()
  res <- gl.filter.rdepth(testset.gl, plot.display = FALSE,
                          plot.dir = tmpdir,
                          plot.file = "chartest_frdepth", verbose = 0)
  saved <- file.path(tmpdir, "chartest_frdepth.RDS")
  expect_true(file.exists(saved))
  unlink(saved)
})

test_that("gl.filter.rdepth returns invisibly", {
  # [approved diff F5] now returns invisibly, matching the other filters.
  v <- withVisible(gl.filter.rdepth(testset.gl, plot.display = FALSE,
                                    verbose = 0))
  expect_false(v$visible)
})

test_that("gl.filter.rdepth progress message wording", {
  # [approved diff F3] message now matches the actual boundary logic.
  out <- capture.output(
    res <- gl.filter.rdepth(testset.gl, lower = 8, upper = 50,
                            plot.display = FALSE, verbose = 2)
  )
  expect_true(any(grepl("rdepth < 8 or > 50", out)))
  expect_false(any(grepl("rdepth <= 8 and >= 50", out)))
})
