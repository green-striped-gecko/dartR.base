# Characterization tests for gl.filter.overshoot
# Baseline snapshotted before review (review-gl.filter.overshoot).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.filter.overshoot removes exactly the overshoot loci", {
  trimmed <- as.character(testset.gl@other$loc.metrics$TrimmedSequence)
  snpos <- testset.gl@other$loc.metrics$SnpPosition
  expected_removed <- sum((snpos + 1) > nchar(trimmed))
  expect_equal(expected_removed, 21)
  res <- gl.filter.overshoot(testset.gl, verbose = 0)
  expect_equal(nLoc(res), nLoc(testset.gl) - expected_removed)
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
  # nothing left that overshoots
  expect_equal(sum((res@other$loc.metrics$SnpPosition + 1) >
                     nchar(as.character(res@other$loc.metrics$TrimmedSequence))),
               0)
})

test_that("gl.filter.overshoot appends one history entry when it filters", {
  h <- length(testset.gl@other$history)
  res <- gl.filter.overshoot(testset.gl, verbose = 0)
  expect_equal(length(res@other$history), h + 1)
})

test_that("gl.filter.overshoot no-op branch returns the input", {
  xno <- testset.gl
  xno@other$loc.metrics$SnpPosition <- rep(0, nLoc(xno))
  out <- capture.output(res <- gl.filter.overshoot(xno, verbose = 0))
  expect_equal(nLoc(res), nLoc(xno))
  # no history appended on a no-op (matches gl.filter.allna, PR #252)
  expect_equal(length(res@other$history), length(xno@other$history))
})

test_that("gl.filter.overshoot leaves the input untouched", {
  x <- testset.gl
  invisible(capture.output(gl.filter.overshoot(x, verbose = 0)))
  expect_identical(x, testset.gl)
})

test_that("gl.filter.overshoot rejects SilicoDArT and missing metrics", {
  expect_error(gl.filter.overshoot(testset.gs, verbose = 0))
  x6 <- testset.gl
  x6@other$loc.metrics$TrimmedSequence <- NULL
  expect_error(gl.filter.overshoot(x6, verbose = 0), "TrimmedSequence")
})

test_that("gl.filter.overshoot no-op message at verbose = 0", {
  # [approved diff F1] no-op message now gates at verbose >= 1.
  xno <- testset.gl
  xno@other$loc.metrics$SnpPosition <- rep(0, nLoc(xno))
  out <- capture.output(res <- gl.filter.overshoot(xno, verbose = 0))
  expect_length(out, 0)
})

test_that("gl.filter.overshoot with NA metric values", {
  # [amended by member direction, 2026-08-31] loci whose overshoot
  # status cannot be assessed (NA TrimmedSequence/SnpPosition) are
  # RETAINED rather than removed (the original review removed them;
  # the coordinator directed retention).
  xna <- testset.gl
  xna@other$loc.metrics$TrimmedSequence[1:5] <- NA
  res <- gl.filter.overshoot(xna, verbose = 0)
  # the 5 NA-metric loci are retained (any genuine overshoots among
  # loci 1:5 cannot be assessed once their metric is NA)
  expect_equal(sum(is.na(res@other$loc.metrics$TrimmedSequence)), 5)
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
  trimmed <- as.character(testset.gl@other$loc.metrics$TrimmedSequence)
  snpos <- testset.gl@other$loc.metrics$SnpPosition + 1
  genuine <- which(snpos > nchar(trimmed))
  expected_removed <- length(setdiff(genuine, 1:5))
  expect_equal(nLoc(res), nLoc(testset.gl) - expected_removed)
})

test_that("gl.filter.overshoot returns invisibly", {
  # [approved diff F2] now returns invisibly, matching the other filters.
  v <- withVisible(gl.filter.overshoot(testset.gl, verbose = 0))
  expect_false(v$visible)
})

test_that("gl.filter.overshoot listing has no stray trailing comma", {
  # [approved diff F4] listing now uses paste(collapse = ", ").
  out <- capture.output(res <- gl.filter.overshoot(testset.gl, verbose = 3))
  listing <- out[grepl("100050384", out)]
  expect_false(grepl(",$", listing))
})
