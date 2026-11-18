# Characterization tests for gl.smearplot
# Baseline snapshotted before review (review-gl.smearplot).
# Assertions marked [approved diff] were flipped in Phase C.

fill_scale <- function(p) {
  Filter(function(s) "fill" %in% s$aesthetics, p$scales$scales)[[1]]
}

test_that("SNP smearplot: visible ggplot, silent at verbose 0, input untouched", {
  pdf(NULL); on.exit(dev.off())
  xcopy <- testset.gl
  o <- capture.output(v <- withVisible(gl.smearplot(testset.gl, verbose = 0)))
  expect_length(o, 0)
  expect_true(v$visible)
  expect_s3_class(v$value, "gg")
  expect_identical(xcopy, testset.gl)
  expect_equal(fill_scale(v$value)$labels,
               c("Homozygote reference", "Heterozygote", "Homozygote alternate"))
})

test_that("plot.display = FALSE renders no page", {
  # [approved diff F1] baseline: plot.display was accepted and guarded
  # but never used - print(p3) was unconditional and a page rendered
  # regardless (file size > 10000).
  tf <- tempfile(fileext = ".pdf")
  pdf(tf)
  invisible(capture.output(
    p <- gl.smearplot(testset.gl, plot.display = FALSE, verbose = 0)))
  dev.off()
  expect_lt(file.info(tf)$size, 10000)   # no page drawn
  unlink(tf)
})

test_that("SilicoDArT legend labels", {
  # [approved diff F2] baseline: labels_silicodart["0"] <- "Absence"
  # appended NAMED elements to an unnamed c("0","1") vector - the legend
  # received c("0","1","Absence","Presence") and displayed "0" and "1".
  pdf(NULL); on.exit(dev.off())
  o <- capture.output(p <- gl.smearplot(testset.gs, verbose = 0))
  lab <- fill_scale(p)$labels
  expect_length(lab, 2)
  expect_equal(unname(lab), c("Absence", "Presence"))
})

test_that("het.only on SilicoDArT", {
  # [approved diff F3] baseline: the het.only palette override happened
  # before the datatype branch, so despite the 'Set to FALSE' warning
  # both states rendered #d3d3d3, and the warning printed at verbose 0.
  pdf(NULL); on.exit(dev.off())
  o <- capture.output(p <- gl.smearplot(testset.gs, het.only = TRUE,
                                        verbose = 0))
  expect_length(o, 0)                                  # warning gated
  expect_equal(fill_scale(p)$palette(2), c("#0000FF", "#FF0000"))
})

test_that("het.only on SNP data uses the het-highlight palette", {
  pdf(NULL); on.exit(dev.off())
  invisible(capture.output(
    p <- gl.smearplot(testset.gl, het.only = TRUE, verbose = 0)))
  expect_equal(fill_scale(p)$palette(3)[2], "#00FFFF")
})

test_that("subset with no homozygote-alternate scores still plots", {
  pdf(NULL); on.exit(dev.off())
  m <- as.matrix(testset.gl)
  no2 <- which(colSums(m == 2, na.rm = TRUE) == 0)
  sub <- testset.gl[, no2[1:5]]
  invisible(capture.output(p <- gl.smearplot(sub, verbose = 0)))
  expect_equal(fill_scale(p)$labels,
               c("Homozygote reference", "Heterozygote"))
})

test_that("loc.order without chromosome info warns and keeps all loci", {
  pdf(NULL); on.exit(dev.off())
  o <- capture.output(p <- gl.smearplot(testset.gl, loc.order = TRUE,
                                        verbose = 2))
  expect_true(any(grepl("chromosome", o)))
  expect_equal(length(unique(p$data$locus)), nLoc(testset.gl))
})

test_that("dendrogram path runs", {
  pdf(NULL); on.exit(dev.off())
  invisible(capture.output(p <- suppressWarnings(
    gl.smearplot(testset.gl[1:15, 1:50], den = TRUE, verbose = 0))))
  expect_false(is.null(p))
})
