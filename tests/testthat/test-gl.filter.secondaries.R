# Characterization tests for gl.filter.secondaries
# Baseline snapshotted before review (review-gl.filter.secondaries).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

clone_of <- function(x) {
  unlist(strsplit(as.character(x@other$loc.metrics$AlleleID),
                  "\\|"))[c(TRUE, FALSE, FALSE)]
}

test_that("gl.filter.secondaries keeps exactly one locus per clone (random)", {
  skip_if_not_installed("dartR.data")
  library(dartR.data)
  set.seed(1)
  res <- gl.filter.secondaries(platypus.gl, method = "random", verbose = 0)
  expect_equal(nLoc(res), 991)
  expect_false(any(duplicated(clone_of(res))))
  expect_setequal(unique(clone_of(platypus.gl)), clone_of(res))
  expect_true(all(ploidy(res) == 2))
  expect_equal(nInd(res), nInd(platypus.gl))
})

test_that("gl.filter.secondaries keeps genotypes and metadata in sync", {
  skip_if_not_installed("dartR.data")
  library(dartR.data)
  set.seed(1)
  res <- gl.filter.secondaries(platypus.gl, method = "random", verbose = 0)
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
  # genotypes of a retained locus are unchanged
  ln <- locNames(res)[1]
  expect_identical(as.matrix(res)[, ln], as.matrix(platypus.gl)[, ln])
})

test_that("gl.filter.secondaries method='best' selection per multi-SNP clone", {
  skip_if_not_installed("dartR.data")
  library(dartR.data)
  res <- gl.filter.secondaries(platypus.gl, method = "best", verbose = 0)
  kept <- as.character(res@other$loc.metrics$AlleleID)
  # [approved diff F1] Pre-fix, selection was lexicographic on the full
  # AlleleID string (RepAvg/AvgPIC sort keys could never engage). Post-fix
  # the highest-RepAvg SNP per clone is retained: clone 45066085 keeps the
  # RepAvg 1.0 SNP (0-25:G>T), clone 45067382 keeps the RepAvg 1.0 SNP
  # (0-64:T>G).
  expect_false("45066085|F|0-23:T>G-23:T>G" %in% kept)
  expect_true("45066085|F|0-25:G>T-25:G>T" %in% kept)
  expect_false("45067382|F|0-10:C>T-10:C>T" %in% kept)
  expect_true("45067382|F|0-64:T>G-64:T>G" %in% kept)
  # every multi-SNP clone retains its best RepAvg (ties broken by AvgPIC)
  b0 <- clone_of(platypus.gl)
  lm0 <- platypus.gl@other$loc.metrics
  for (cl in names(table(b0))[table(b0) > 1]) {
    rows <- lm0[b0 == cl, c("AlleleID", "RepAvg", "AvgPIC")]
    kept_id <- kept[kept %in% as.character(rows$AlleleID)]
    kept_rep <- rows$RepAvg[as.character(rows$AlleleID) == kept_id]
    expect_equal(kept_rep, max(rows$RepAvg))
  }
})

test_that("gl.filter.secondaries locus order in the output", {
  skip_if_not_installed("dartR.data")
  library(dartR.data)
  # [approved diff F2] Pre-fix the output locus order was permuted for both
  # methods, even when no secondaries existed. Post-fix the input locus
  # order is preserved; a no-op filter returns the loci unchanged.
  set.seed(1)
  res_ts <- gl.filter.secondaries(testset.gl, verbose = 0)
  expect_equal(nLoc(res_ts), 255)
  expect_identical(locNames(res_ts), locNames(testset.gl))
  set.seed(1)
  res_pl <- gl.filter.secondaries(platypus.gl, verbose = 0)
  in_order <- locNames(platypus.gl)
  expect_identical(locNames(res_pl), in_order[in_order %in% locNames(res_pl)])
})

test_that("gl.filter.secondaries verbosity and warnings", {
  skip_if_not_installed("dartR.data")
  library(dartR.data)
  out <- capture.output(res <- gl.filter.secondaries(platypus.gl, verbose = 0))
  expect_length(out, 0)
  # [approved diff F3] invalid-method warning now gates at verbose >= 1
  out2 <- capture.output(
    res2 <- gl.filter.secondaries(platypus.gl, method = "bogus", verbose = 0)
  )
  expect_length(out2, 0)
  # invalid method behaves as 'random'
  expect_equal(nLoc(res2), 991)
})

test_that("gl.filter.secondaries appends one history entry and leaves input untouched", {
  skip_if_not_installed("dartR.data")
  library(dartR.data)
  x <- platypus.gl
  h <- length(x@other$history)
  set.seed(1)
  res <- gl.filter.secondaries(x, verbose = 0)
  expect_identical(x, platypus.gl)
  expect_equal(length(res@other$history), h + 1)
})

test_that("gl.filter.secondaries returns invisibly", {
  skip_if_not_installed("dartR.data")
  library(dartR.data)
  # [approved diff F4] now returns invisibly, matching the other filters.
  set.seed(1)
  v <- withVisible(gl.filter.secondaries(platypus.gl, verbose = 0))
  expect_false(v$visible)
})

test_that("gl.filter.secondaries rejects SilicoDArT data", {
  expect_error(gl.filter.secondaries(testset.gs, verbose = 0))
})
