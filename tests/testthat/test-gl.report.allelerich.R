# Characterization tests for gl.report.allelerich
# Baseline snapshotted before review (review-gl.report.allelerich).
# Assertions marked [approved diff] were flipped in Phase C.

gl4 <- gl.filter.allna(gl.keep.pop(testset.gl,
        pop.list = c("EmmacBrisWive", "EmmacBurdMist",
                     "EmmacClarJack", "EmmacRussEube"),
        verbose = 0), verbose = 0)

test_that("richness values match an independent recomputation", {
  pdf(NULL); on.exit(dev.off())
  invisible(capture.output(res <- gl.report.allelerich(gl4,
        plot.display = FALSE, verbose = 0)))
  expect_named(res, c("Allelic Richness per site",
                      "Allelic Richness per population",
                      "Raw reference allele count",
                      "Raw alternate allele count"))
  fn <- res[["Allelic Richness per population"]]
  means <- fn$mean_corrected_richness
  names(means) <- fn$pop
  expect_equal(round(means[["EmmacBrisWive"]], 5), 1.01484)
  expect_equal(round(means[["EmmacBurdMist"]], 5), 1.00990)
  expect_equal(round(means[["EmmacClarJack"]], 5), 1.01150)
  expect_equal(round(means[["EmmacRussEube"]], 5), 1.02391)
  expect_equal(fn$popsize, c(10, 10, 5, 10))
  expect_true("SD" %in% colnames(fn))          # default error.bar
})

test_that("silence and visibility at verbose 0", {
  pdf(NULL); on.exit(dev.off())
  xcopy <- gl4
  o <- capture.output(v <- withVisible(
    gl.report.allelerich(gl4, plot.display = FALSE, verbose = 0)))
  # [approved diff F1] baseline: the lazy default gl.colors("dis")
  # printed the 3-line gl.colors banner at every verbosity.
  expect_length(o, 0)
  expect_false(v$visible)
  expect_identical(xcopy, gl4)
})

test_that("plot rendering at verbose 0", {
  # BASELINE: no verbose==0 -> plot.display <- FALSE guard; the plot
  # renders at verbose 0.
  t0 <- tempfile(fileext = ".pdf"); pdf(t0); dev.off()
  empty <- file.info(t0)$size
  tf <- tempfile(fileext = ".pdf"); pdf(tf)
  invisible(capture.output(g <- gl.report.allelerich(gl4, verbose = 0)))
  dev.off()
  # [approved diff F1] baseline: a page was drawn at verbose 0.
  expect_lt(file.info(tf)$size, empty + 500)
  unlink(c(t0, tf))
})

test_that("error.bar validation", {
  pdf(NULL); on.exit(dev.off())
  # [approved diff F2] baseline: an unknown error.bar crashed
  # downstream with "object 'max_val' not found". Now coerced to SD.
  invisible(capture.output(res <- gl.report.allelerich(gl4,
        error.bar = "bogus", plot.display = FALSE, verbose = 0)))
  expect_true("SD" %in% colnames(res[["Allelic Richness per population"]]))
})

test_that("bootstrap path adds CI columns", {
  pdf(NULL); on.exit(dev.off())
  skip_if_not_installed("boot")
  skip_if_not_installed("Rcpp")
  o <- capture.output(res <- suppressWarnings(
    gl.report.allelerich(gl4[, 1:60], nboots = 15, CI.type = "perc",
                         error.bar = "SD", plot.display = FALSE,
                         verbose = 0)))
  cn <- colnames(res[["Allelic Richness per population"]])
  expect_true(all(c("LCI", "HCI") %in% cn))
  # BASELINE: the user's error.bar = "SD" is silently overridden to CI.
  expect_false("SD" %in% cn)
})
