# Characterization tests for gl.filter.hamming
# Baseline snapshotted before review (dev at 984ab05). Expectations marked
# CURRENT DEFECT pin today's behaviour; they flip only if an approved finding
# changes it.

test_that("default filter removes 11 of 1000 platypus loci, metadata in sync", {
  skip_if_not_installed("Rcpp")
  x <- platypus.gl
  f <- gl.filter.hamming(x, verbose = 0)
  expect_equal(nLoc(f), 989L)
  expect_equal(nrow(f@other$loc.metrics), 989L)
  i <- match(locNames(f), locNames(x))
  expect_identical(f@other$loc.metrics, x@other$loc.metrics[i, , drop = FALSE])
  expect_equal(length(f@other$history), length(x@other$history) + 1L)
  expect_equal(unique(ploidy(f)), 2L)
})

test_that("no two kept comparable loci are within threshold (property)", {
  skip_if_not_installed("Rcpp")
  f <- gl.filter.hamming(platypus.gl, threshold = 3, verbose = 0)
  s <- substr(toupper(as.character(f@other$loc.metrics$TrimmedSequence)),
              6, 55)
  s <- s[nchar(s) == 50]
  m <- do.call(rbind, strsplit(s, ""))
  md <- Inf
  for (i in 1:(nrow(m) - 1)) {
    md <- min(md, min(colSums(t(m[(i + 1):nrow(m), , drop = FALSE]) != m[i, ])))
  }
  expect_gt(md, 3)
})

test_that("removal counts across thresholds", {
  skip_if_not_installed("Rcpp")
  x <- platypus.gl
  n <- function(...) nLoc(x) - nLoc(gl.filter.hamming(x, ..., verbose = 0))
  expect_equal(n(threshold = 0), 8L)
  expect_equal(n(threshold = 3), 11L)
  expect_equal(n(threshold = 10), 20L)
})

test_that("SilicoDArT input keeps ploidy 1", {
  skip_if_not_installed("Rcpp")
  g <- gl.filter.hamming(testset.gs, verbose = 0)
  expect_equal(nLoc(g), 237L)
  expect_equal(unique(ploidy(g)), 1L)
})

test_that("proportion-style threshold errors", {
  skip_if_not_installed("Rcpp")
  expect_error(gl.filter.hamming(platypus.gl, threshold = 0.2, verbose = 0),
               "not a proportion")
})

test_that("invalid arguments stop with an error", {
  skip_if_not_installed("Rcpp")
  x <- platypus.gl
  f <- function(...) gl.filter.hamming(x, ..., verbose = 0)
  # Change 2 (F2): threshold must be one whole number >= 0; 2.9 ran as 2
  expect_error(f(threshold = 2.9), "whole number")
  expect_error(f(threshold = NA), "whole number")
  expect_error(f(threshold = c(1, 3)), "whole number")
  expect_error(f(threshold = -1), "whole number")
  # Change 1 (F1): threshold >= min.length removed 948 of 1000 loci
  expect_error(f(threshold = 50), "smaller than")
  # Change 3 (F3): invalid rs / min.length compared nothing, silently
  expect_error(f(rs = -3), "rs")
  expect_error(f(min.length = 0), "min.length")
  expect_error(f(min.length = 2.5), "min.length")
})

test_that("nothing comparable warns at verbose 1 and returns x", {
  skip_if_not_installed("Rcpp")
  x <- platypus.gl
  out <- capture.output(r <- gl.filter.hamming(x, min.length = 100,
                                              verbose = 1))
  expect_true(any(grepl("fewer than two loci", out)))
  expect_equal(nLoc(r), nLoc(x))
})
