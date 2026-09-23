# Characterization tests for gl.report.hamming
# Baseline snapshotted before review (dev at 984ab05). Expectations marked
# CURRENT DEFECT pin today's behaviour; they flip only if an approved finding
# changes it.

test_that("report returns the loci-removed table", {
  skip_if_not_installed("Rcpp")
  x <- platypus.gl
  # Change 4 (F4): previously returned x unchanged
  r <- gl.report.hamming(x, verbose = 0)
  expect_s3_class(r, "data.frame")
  expect_named(r, c("Threshold", "Removed", "Percent.removed", "Retained",
                    "Percent.retained"))
  expect_equal(r$Threshold, 0:10)
  expect_equal(r$Removed, c(8, 9, 10, 11, 12, 12, 15, 15, 17, 19, 20))
  expect_equal(r$Retained, nLoc(x) - r$Removed)
})

test_that("loci-removed table matches gl.filter.hamming exactly", {
  skip_if_not_installed("Rcpp")
  x <- platypus.gl
  out <- capture.output(
    gl.report.hamming(x, verbose = 2, plot.display = FALSE)
  )
  rows <- grep("^ +[0-9]+ +[0-9]+ ", out, value = TRUE)
  tab <- read.table(text = rows)
  expect_equal(tab$V1, 0:10)
  expect_equal(tab$V2, c(8, 9, 10, 11, 12, 12, 15, 15, 17, 19, 20))
  for (k in c(0, 3, 10)) {
    f <- gl.filter.hamming(x, threshold = k, verbose = 0)
    expect_equal(tab$V2[tab$V1 == k], nLoc(x) - nLoc(f))
  }
})

test_that("verbose 1 prints only start and end; table is returned", {
  skip_if_not_installed("Rcpp")
  out <- capture.output(r <- gl.report.hamming(platypus.gl, verbose = 1))
  expect_length(out, 2L)
  expect_equal(nrow(r), 11L)
})

test_that("argument checks", {
  skip_if_not_installed("Rcpp")
  x <- platypus.gl
  expect_error(gl.report.hamming(x, threshold = 0.2, verbose = 0),
               "not a proportion")
  # Change 3 (F3): shared validation message
  expect_error(gl.report.hamming(x, rs = -3, verbose = 0), "rs")
  # Changes 1, 2 (F1, F2)
  expect_error(gl.report.hamming(x, threshold = 2.9, verbose = 0),
               "whole number")
  expect_error(gl.report.hamming(x, threshold = 50, verbose = 0),
               "smaller than")
})
