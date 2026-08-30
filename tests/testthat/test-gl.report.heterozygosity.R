# Characterization tests for gl.report.heterozygosity
# Baseline snapshotted before review (review-gl.report.heterozygosity).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("gl.report.heterozygosity pop statistics match hand computation", {
  pdf(NULL); on.exit(dev.off())
  r <- as.data.frame(gl.report.heterozygosity(testset.gl,
                                              plot.display = FALSE,
                                              verbose = 0))
  p1 <- "EmmacBurdMist"
  m <- as.matrix(testset.gl[as.character(pop(testset.gl)) == p1, ])
  Ho_hand <- mean(colMeans(m == 1, na.rm = TRUE), na.rm = TRUE)
  n0 <- colSums(m == 0, na.rm = TRUE)
  n1 <- colSums(m == 1, na.rm = TRUE)
  n2 <- colSums(m == 2, na.rm = TRUE)
  p <- (2 * n0 + n1) / (2 * (n0 + n1 + n2))
  He_hand <- mean(1 - (p^2 + (1 - p)^2), na.rm = TRUE)
  row <- r[r$pop == p1, ]
  expect_equal(row$Ho, round(Ho_hand, 6), tolerance = 1e-4)
  expect_equal(row$He, round(He_hand, 6), tolerance = 1e-3)
  expect_true(all(c("uHe", "FIS", "polyLoc", "monoLoc", "all_NALoc")
                  %in% names(r)))
})

test_that("gl.report.heterozygosity is silent at verbose 0 and invisible", {
  pdf(NULL); on.exit(dev.off())
  o1 <- capture.output(
    r1 <- gl.report.heterozygosity(testset.gl, plot.display = FALSE,
                                   verbose = 0))
  o2 <- capture.output(
    r2 <- gl.report.heterozygosity(testset.gl, method = "ind",
                                   plot.display = FALSE, verbose = 0))
  expect_length(o1, 0)
  expect_length(o2, 0)
  v <- withVisible(gl.report.heterozygosity(testset.gl,
                                            plot.display = FALSE,
                                            verbose = 0))
  expect_false(v$visible)
  xcopy <- testset.gl
  invisible(capture.output(gl.report.heterozygosity(xcopy,
                                                    plot.display = FALSE,
                                                    verbose = 0)))
  expect_identical(xcopy, testset.gl)
})

test_that("method='ind' with display off at verbose 3", {
  # BASELINE (pre-fix): crashes with "object 'outliers' not found" — the
  # outlier table is built only inside the plotting block.
  pdf(NULL); on.exit(dev.off())
  expect_error(
    capture.output(gl.report.heterozygosity(testset.gl, method = "ind",
                                            plot.display = FALSE,
                                            verbose = 3)),
    "outliers"
  )
})

test_that("subsample.pop with method='ind'", {
  # BASELINE (pre-fix): crashes with "object 'res_sub' not found" — the
  # subsample results exist only under method='pop'.
  pdf(NULL); on.exit(dev.off())
  expect_error(
    capture.output(gl.report.heterozygosity(testset.gl, method = "ind",
                                            subsample.pop = TRUE,
                                            plot.display = FALSE,
                                            verbose = 0)),
    "res_sub"
  )
})

test_that("subsample.pop crashes when any population is below n.limit", {
  # BASELINE (pre-fix): utils.subsample.pop stores NA for below-limit
  # populations and rbindlist() rejects it — testset.gl (many pops < 10)
  # crashes although n.limit is documented as a skip threshold.
  pdf(NULL); on.exit(dev.off())
  expect_error(
    capture.output(gl.report.heterozygosity(testset.gl,
                                            subsample.pop = TRUE,
                                            plot.display = FALSE,
                                            verbose = 0)),
    "plain list"
  )
})

test_that("subsample.pop works when all populations qualify (platypus)", {
  pdf(NULL); on.exit(dev.off())
  res <- gl.report.heterozygosity(platypus.gl, subsample.pop = TRUE,
                                  plot.display = FALSE, verbose = 0)
  # BASELINE (pre-fix): an unnamed 2-element list
  expect_true(is.list(res) && !is.data.frame(res))
  expect_length(res, 2)
  expect_null(names(res))
})

test_that("warnings print at verbose 0", {
  # BASELINE (pre-fix): the method-coercion and negative-n.invariant
  # warnings are ungated.
  pdf(NULL); on.exit(dev.off())
  o5 <- capture.output(
    r5 <- gl.report.heterozygosity(testset.gl, method = "bogus",
                                   plot.display = FALSE, verbose = 0))
  expect_gt(length(o5), 0)
  o6 <- capture.output(
    r6 <- gl.report.heterozygosity(testset.gl, n.invariant = -1,
                                   plot.display = FALSE, verbose = 0))
  expect_gt(length(o6), 0)
})
