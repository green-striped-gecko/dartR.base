# Characterization tests for gl.filter.taglength
# Baseline snapshotted before review (review-gl.filter.taglength).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.filter.taglength retains loci within [lower, upper] inclusive", {
  ncht <- nchar(as.character(testset.gl@other$loc.metrics$TrimmedSequence))
  expected <- sum(ncht >= 60 & ncht <= 69)
  res <- gl.filter.taglength(testset.gl, lower = 60, verbose = 0)
  expect_equal(nLoc(res), expected)
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
  out.len <- nchar(as.character(res@other$loc.metrics$TrimmedSequence))
  expect_true(all(out.len >= 60 & out.len <= 69))
  # boundary loci retained
  expect_equal(sum(out.len == 60), sum(ncht == 60))
})

test_that("gl.filter.taglength works on SilicoDArT with TrimmedSequence", {
  nchs <- nchar(as.character(testset.gs@other$loc.metrics$TrimmedSequence))
  res <- gl.filter.taglength(testset.gs, lower = 60, verbose = 0)
  expect_equal(nLoc(res), sum(nchs >= 60 & nchs <= 69))
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
})

test_that("gl.filter.taglength appends history and leaves the input untouched", {
  x <- testset.gl
  h <- length(x@other$history)
  res <- gl.filter.taglength(x, verbose = 0)
  expect_identical(x, testset.gl)
  expect_equal(length(res@other$history), h + 1)
})

test_that("gl.filter.taglength errors when TrimmedSequence is absent", {
  x <- testset.gl
  x@other$loc.metrics$TrimmedSequence <- NULL
  expect_error(gl.filter.taglength(x, verbose = 0), "TrimmedSequence")
})

test_that("gl.filter.taglength swaps reversed thresholds", {
  out <- capture.output(
    res <- gl.filter.taglength(testset.gl, lower = 69, upper = 60, verbose = 0)
  )
  ncht <- nchar(as.character(testset.gl@other$loc.metrics$TrimmedSequence))
  expect_equal(nLoc(res), sum(ncht >= 60 & ncht <= 69))
})

test_that("gl.filter.taglength with NA tag lengths", {
  # [approved diff F1] Pre-fix, NA TrimmedSequence values silently DESYNCED
  # genotypes and metadata (185 loci vs 195 metric rows). Post-fix NA-length
  # loci are removed cleanly and sync holds.
  xna <- testset.gl
  xna@other$loc.metrics$TrimmedSequence[1:10] <- NA
  res <- gl.filter.taglength(xna, lower = 60, verbose = 0)
  expect_equal(nLoc(res), 185)
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
  expect_false(any(is.na(res@other$loc.metrics$TrimmedSequence)))
})

test_that("gl.filter.taglength warning gating at verbose = 0", {
  # [approved diff F3] warnings now gate at verbose >= 1.
  out <- capture.output(
    res <- gl.filter.taglength(testset.gl, lower = 69, upper = 60, verbose = 0)
  )
  expect_length(out, 0)
})

test_that("gl.filter.taglength returns invisibly", {
  # [approved diff F4] now returns invisibly, matching the other filters.
  v <- withVisible(gl.filter.taglength(testset.gl, verbose = 0))
  expect_false(v$visible)
})

test_that("gl.filter.taglength progress message wording", {
  # [approved diff F2] message now matches the actual behaviour.
  out <- capture.output(
    res <- gl.filter.taglength(testset.gl, lower = 60, verbose = 2)
  )
  expect_true(any(grepl("tag length < 60 or > 69", out)))
  expect_false(any(grepl("taglength between", out)))
})
