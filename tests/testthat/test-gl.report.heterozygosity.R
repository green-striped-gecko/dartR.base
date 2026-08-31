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
  # [approved diff F2] baseline: crashed with "object 'outliers' not
  # found" — the outlier table was built only inside the plotting block.
  pdf(NULL); on.exit(dev.off())
  o <- capture.output(r <- gl.report.heterozygosity(testset.gl,
        method = "ind", plot.display = FALSE, verbose = 3))
  expect_s3_class(r, "data.frame")
  expect_true(any(grepl("[Oo]utliers", o)))
})

test_that("subsample.pop with method='ind'", {
  # [approved diff F3] baseline: crashed with "object 'res_sub' not
  # found". Now a gated warning; subsample.pop is ignored and the plain
  # dataframe is returned.
  pdf(NULL); on.exit(dev.off())
  o <- capture.output(r <- gl.report.heterozygosity(testset.gl,
        method = "ind", subsample.pop = TRUE, plot.display = FALSE,
        verbose = 0))
  expect_length(o, 0)
  expect_s3_class(r, "data.frame")
  o1 <- capture.output(r1 <- gl.report.heterozygosity(testset.gl,
        method = "ind", subsample.pop = TRUE, plot.display = FALSE,
        verbose = 1))
  expect_true(any(grepl("ignored", o1)))
})

test_that("subsample.pop skips populations below n.limit", {
  # [approved diff F1] baseline: utils.subsample.pop stored NA for
  # below-limit populations and rbindlist() rejected it — testset.gl
  # (many pops < 10) crashed although n.limit is documented as a skip
  # threshold.
  pdf(NULL); on.exit(dev.off())
  res <- gl.report.heterozygosity(testset.gl, subsample.pop = TRUE,
                                  plot.display = FALSE, verbose = 0)
  expect_named(res, c("subsample", "results"))
  qualifying <- names(which(table(pop(testset.gl)) >= 10))
  expect_setequal(unique(res$subsample$pop), qualifying)
  expect_equal(nrow(res$results), nPop(testset.gl))
})

test_that("subsample.pop works when all populations qualify (platypus)", {
  pdf(NULL); on.exit(dev.off())
  res <- gl.report.heterozygosity(platypus.gl, subsample.pop = TRUE,
                                  plot.display = FALSE, verbose = 0)
  # [approved diff F5] baseline: an unnamed 2-element list.
  expect_true(is.list(res) && !is.data.frame(res))
  expect_length(res, 2)
  expect_named(res, c("subsample", "results"))
  expect_setequal(unique(res$subsample$pop), popNames(platypus.gl))
})

test_that("warnings gated at verbose 0", {
  # [approved diff F4] baseline: the method-coercion and
  # negative-n.invariant warnings printed at verbose 0.
  pdf(NULL); on.exit(dev.off())
  o5 <- capture.output(
    r5 <- gl.report.heterozygosity(testset.gl, method = "bogus",
                                   plot.display = FALSE, verbose = 0))
  expect_length(o5, 0)
  o6 <- capture.output(
    r6 <- gl.report.heterozygosity(testset.gl, n.invariant = -1,
                                   plot.display = FALSE, verbose = 0))
  expect_length(o6, 0)
})
