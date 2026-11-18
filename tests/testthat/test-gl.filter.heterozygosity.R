# Characterization tests for gl.filter.heterozygosity
# Baseline snapshotted before review (review-gl.filter.heterozygosity).
# Assertions marked [approved diff] were flipped in Phase C.

make_allna <- function() {
  m2 <- matrix(c(0,1,2,0, NA,NA,NA,NA, 2,1,0,1), nrow = 3, byrow = TRUE,
               dimnames = list(c("IND_A","IND_ALLNA","IND_C"),
                               paste0("L",1:4)))
  gx <- new("genlight", gen = m2, ploidy = 2)
  pop(gx) <- factor(rep("p1", 3))
  gl.compliance.check(gx, verbose = 0)
}

test_that("gl.filter.heterozygosity retains the right individuals", {
  m <- as.matrix(testset.gl)
  c.na <- rowSums(is.na(m))
  c.hets <- rowSums(m == 1, na.rm = TRUE) / (nLoc(testset.gl) - c.na)
  expected <- sum(c.hets >= 0 & c.hets <= 0.06)
  out <- capture.output(
    res <- gl.filter.heterozygosity(testset.gl, t.upper = 0.06, verbose = 0)
  )
  expect_equal(nInd(res), expected)
  expect_equal(nLoc(res), nLoc(testset.gl))
  expect_equal(nrow(res@other$ind.metrics), nInd(res))
  expect_identical(as.character(res@other$ind.metrics$id), indNames(res))
})

test_that("gl.filter.heterozygosity history +1, input untouched", {
  h <- length(testset.gl@other$history)
  out <- capture.output(
    res <- gl.filter.heterozygosity(testset.gl, t.upper = 0.06, verbose = 0)
  )
  expect_equal(length(res@other$history), h + 1)
  xcopy <- testset.gl
  invisible(capture.output(
    gl.filter.heterozygosity(xcopy, t.upper = 0.06, verbose = 0)))
  expect_identical(xcopy, testset.gl)
})

test_that("locus-metric flags after individual removal", {
  # [approved diff F1] flags are now reset when individuals are removed,
  # matching the gl.drop.ind precedent, at every verbosity.
  res <- gl.filter.heterozygosity(testset.gl, t.upper = 0.06, verbose = 0)
  expect_gt(nInd(testset.gl) - nInd(res), 0)  # something was removed
  expect_false(res@other$loc.metrics.flags$CallRate)
  expect_false(res@other$loc.metrics.flags$FreqHets)
  # no removal -> flags untouched
  res2 <- gl.filter.heterozygosity(testset.gl, t.upper = 1, verbose = 0)
  expect_equal(nInd(res2), nInd(testset.gl))
  expect_true(res2@other$loc.metrics.flags$CallRate)
})

test_that("all-NA individuals", {
  # [approved diff F2] individuals whose heterozygosity cannot be computed
  # are removed cleanly instead of crashing.
  gx <- make_allna()
  res <- gl.filter.heterozygosity(gx, t.upper = 0.7, verbose = 0)
  expect_equal(nInd(res), 2)
  expect_false("IND_ALLNA" %in% indNames(res))
})

test_that("monomorphs warning at verbose 0", {
  # [approved diff F3] the warning now gates at verbose >= 2.
  out <- capture.output(
    res <- gl.filter.heterozygosity(testset.gl, t.upper = 0.06, verbose = 0)
  )
  expect_length(out, 0)
})

test_that("flag-less object", {
  # [approved diff F3] isFALSE() guard: flag-less objects no longer crash.
  gx <- make_allna()
  gx@other$loc.metrics.flags <- NULL
  res <- gl.filter.heterozygosity(gx, verbose = 0)
  expect_equal(nInd(res), 2)
})

test_that("t.lower > t.upper", {
  # [approved diff F4] reversed thresholds are swapped with a gated
  # warning (silent at verbose 0), matching gl.filter.taglength.
  out <- capture.output(
    res <- gl.filter.heterozygosity(testset.gl, t.lower = 0.5,
                                    t.upper = 0.01, verbose = 0)
  )
  expect_length(out, 0)
  m <- as.matrix(testset.gl)
  c.na <- rowSums(is.na(m))
  c.hets <- rowSums(m == 1, na.rm = TRUE) / (nLoc(testset.gl) - c.na)
  expect_equal(nInd(res), sum(c.hets >= 0.01 & c.hets <= 0.5))
})

test_that("t.lower error message text", {
  # [approved diff F4] the t.lower range error now names t.lower.
  err <- tryCatch(
    gl.filter.heterozygosity(testset.gl, t.lower = -0.5, verbose = 0),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("t.lower", err))
})

test_that("returns invisibly", {
  # [approved diff F5] now returns invisibly.
  v <- withVisible(gl.filter.heterozygosity(testset.gl, t.upper = 0.06,
                                            verbose = 0))
  expect_false(v$visible)
})
