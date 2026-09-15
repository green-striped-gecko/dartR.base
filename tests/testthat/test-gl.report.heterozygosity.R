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

# ---- PR #229 (Carlo Pacioni) re-applied on the reviewed code ----------------

test_that("uHe and FIS use Nei's per-locus sample-size correction", {
  pdf(NULL); on.exit(dev.off())
  r <- as.data.frame(gl.report.heterozygosity(testset.gl,
                                              plot.display = FALSE,
                                              verbose = 0))
  p1 <- "EmmacBurdMist"
  m <- as.matrix(testset.gl[as.character(pop(testset.gl)) == p1, ])
  m <- m[, colSums(!is.na(m)) > 0, drop = FALSE]      # all-NA loci excluded
  n_l <- colSums(!is.na(m))                           # genotyped per locus
  q <- colMeans(m, na.rm = TRUE) / 2
  He_l <- 2 * q * (1 - q)
  uHe_l <- (2 * n_l / (2 * n_l - 1)) * He_l
  Ho_l <- colMeans(m == 1, na.rm = TRUE)
  row <- r[r$pop == p1, ]
  expect_equal(row$uHe, mean(uHe_l), tolerance = 1e-5)
  expect_equal(row$FIS, mean(1 - Ho_l / uHe_l, na.rm = TRUE),
               tolerance = 1e-5)
})

test_that("adjusted SEs use scored loci + n.invariant as the sample size", {
  pdf(NULL); on.exit(dev.off())
  r <- as.data.frame(gl.report.heterozygosity(testset.gl, n.invariant = 100,
                                              plot.display = FALSE,
                                              verbose = 0))
  expect_equal(r$Ho.adjSE, r$Ho.adjSD / sqrt(r$n.Loc + 100), tolerance = 1e-4)
  expect_equal(r$He.adjSE, r$He.adjSD / sqrt(r$n.Loc + 100), tolerance = 1e-4)
})

test_that("the bootstrap statistic is the reported point estimate", {
  pdf(NULL); on.exit(dev.off())
  r <- as.data.frame(gl.report.heterozygosity(testset.gl,
                                              plot.display = FALSE,
                                              verbose = 0))
  p1 <- "EmmacBurdMist"
  m <- as.matrix(testset.gl[as.character(pop(testset.gl)) == p1, ])
  stat <- pop.het(m, indices = seq_len(nrow(m)), n.invariant = 0,
                  boot_method = "ind", aHet = FALSE)
  row <- r[r$pop == p1, ]
  expect_equal(unname(stat), c(row$Ho, row$He, row$uHe, row$FIS),
               tolerance = 1e-5)
})

test_that("single-individual populations survive the plain and bootstrap paths", {
  # baseline: 'x' must be an array of at least two dimensions (the 1-row
  # matrix collapsed to a vector) in both paths; testset.gl has two such
  # populations
  pdf(NULL); on.exit(dev.off())
  one <- names(which(table(pop(testset.gl)) == 1))[1]
  x <- gl.keep.pop(testset.gl, pop.list = c(one, "EmmacBurdMist"),
                   verbose = 0)
  r <- gl.report.heterozygosity(x, plot.display = FALSE, verbose = 0)
  expect_s3_class(r, "data.frame")
  expect_equal(r$n.Ind[r$pop == one], 1)
  set.seed(1)
  suppressWarnings(
    rb <- gl.report.heterozygosity(x, nboots = 30, CI.type = "perc",
                                   plot.display = FALSE, verbose = 0))
  expect_s3_class(rb, "data.frame")
  expect_true(all(c("HoLCI", "HoHCI") %in% names(rb)))
  expect_false(anyNA(rb[rb$pop == "EmmacBurdMist", c("HoLCI", "HoHCI")]))
})
