# Characterization tests for gl.report.basics
# Baseline snapshotted before review (review-gl.report.basics).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("silent at verbose 0; invisible NULL; input untouched", {
  xcopy <- testset.gl
  o <- capture.output(v <- withVisible(gl.report.basics(testset.gl,
                                                        verbose = 0)))
  expect_length(o, 0)
  expect_false(v$visible)
  expect_null(v$value)
  expect_identical(xcopy, testset.gl)
})

test_that("SNP report at verbose 1 carries the key statistics", {
  o <- capture.output(gl.report.basics(testset.gl, verbose = 1))
  expect_true(any(grepl("Loci: +255", o)))
  expect_true(any(grepl("Individuals: +250", o)))
  expect_true(any(grepl("Populations: +30", o)))
  # consistent with gl.filter.monomorphs (144 on testset.gl)
  expect_true(any(grepl("Monomorphic Loci: +144", o)))
  expect_true(any(grepl("Loci all NA: +0", o)))
})

test_that("SilicoDArT data", {
  # [approved diff F1] baseline: crashed - the composition table
  # hard-coded four column names onto a three-class table.
  o <- capture.output(gl.report.basics(testset.gs, verbose = 1))
  expect_true(any(grepl("SilicoDArT", o)))
  expect_true(any(grepl("Monomorphic Loci:", o)))
})

test_that("SNP subset with a genotype class absent", {
  # [approved diff F1] baseline: crashed for the same reason.
  m <- as.matrix(testset.gl)
  no2 <- which(colSums(m == 2, na.rm = TRUE) == 0)
  sub <- testset.gl[, no2[1:10]]
  o <- capture.output(gl.report.basics(sub, verbose = 1))
  expect_true(any(grepl("Loci: +10", o)))
})

make_gg <- function() {
  set.seed(11)
  mm <- matrix(sample(c(0, 1, 2), 20 * 30, replace = TRUE), nrow = 20)
  mm[3, ] <- NA
  gg <- new("genlight", gen = mm, ploidy = 2)
  indNames(gg) <- paste0("ind", 1:20)
  locNames(gg) <- paste0("loc", 1:30)
  pop(gg) <- factor(rep(c("A", "B"), each = 10))
  suppressWarnings(gl.compliance.check(gg, verbose = 0))
}

test_that("all-NA individuals listing", {
  gg <- make_gg()
  o <- suppressWarnings(capture.output(gl.report.basics(gg, verbose = 1)))
  expect_true(any(grepl("Individuals all NA: +1", o)))
  iline <- grep("Individuals all NA", o)[1]
  # [approved diff F2] baseline: the listing printed the whole
  # NA-padded array ("NA NA ind3 NA ...").
  expect_true(grepl("^ind3", trimws(o[iline + 1])))
  expect_false(grepl("NA NA", o[iline + 1]))
})

test_that("missing rdepth metric", {
  gg <- make_gg()
  gg@other$loc.metrics$rdepth <- NULL
  # [approved diff F3] baseline: raw mean.default warning and NA.
  expect_no_warning(o <- capture.output(gl.report.basics(gg, verbose = 1)))
  expect_true(any(grepl("not available", o)))
})
