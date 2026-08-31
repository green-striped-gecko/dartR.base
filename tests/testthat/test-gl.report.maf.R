# Characterization tests for gl.report.maf
# Baseline snapshotted before review (review-gl.report.maf).
# Assertions marked [approved diff] were flipped in Phase C.

gl4 <- gl.keep.pop(testset.gl,
                   pop.list = c("EmmacBrisWive", "EmmacBurdMist",
                                "EmmacClarJack", "EmmacRussEube"),
                   verbose = 0)

test_that("MAF values match a hand computation; return structure", {
  pdf(NULL); on.exit(dev.off())
  invisible(capture.output(res <- gl.report.maf(gl4, plot.display = FALSE,
                                                verbose = 0)))
  expect_s3_class(res, "data.frame")
  expect_equal(dim(res), c(255, 4))
  expect_setequal(colnames(res),
                  c("EmmacBrisWive", "EmmacBurdMist",
                    "EmmacClarJack", "EmmacRussEube"))
  m <- as.matrix(gl4[as.character(pop(gl4)) == "EmmacBurdMist", ])
  q <- colMeans(m, na.rm = TRUE) / 2
  mymaf <- pmin(q, 1 - q)
  fn <- res[["EmmacBurdMist"]]
  keep <- !is.na(mymaf) & !is.na(fn)
  expect_true(all(abs(mymaf[keep] - fn[keep]) < 1e-6))
})

test_that("verbose 0 silence and visibility", {
  pdf(NULL); on.exit(dev.off())
  xcopy <- gl4
  o <- capture.output(v <- withVisible(gl.report.maf(gl4,
        plot.display = FALSE, verbose = 0)))
  # [approved diff F1] baseline: 77 lines printed at verbose 0
  # (per-population stats via a hardcoded verbose = 3 inside the
  # per-pop function, ungated overall stats, ungated quantile table).
  expect_length(o, 0)
  expect_false(v$visible)
  expect_identical(xcopy, gl4)
})

test_that("single Starting banner", {
  pdf(NULL); on.exit(dev.off())
  o <- capture.output(r <- gl.report.maf(gl4, plot.display = FALSE,
                                         verbose = 2))
  # [approved diff F2] baseline: a second handwritten FLAG SCRIPT START
  # block duplicated the banner (reading a package-level 'build'
  # constant by lexical accident).
  expect_equal(sum(grepl("Starting", o)), 1)
})

test_that("as.pop reassignment works and restores", {
  pdf(NULL); on.exit(dev.off())
  sub <- testset.gl[1:40, ]
  invisible(capture.output(r <- gl.report.maf(sub, as.pop = "sex",
        plot.display = FALSE, verbose = 0)))
  expect_setequal(colnames(r), c("Female", "Male", "Unknown"))
  expect_error(capture.output(gl.report.maf(gl4, as.pop = "bogus",
        plot.display = FALSE, verbose = 0)))
})

test_that("singleton populations dropped from the analysis", {
  pdf(NULL); on.exit(dev.off())
  # [approved diff F1+F3] baseline: the message printed at verbose 0
  # and reported the pop INDEX (e.g. "23"), not the population name.
  o0 <- capture.output(r <- gl.report.maf(testset.gl, plot.display = FALSE,
                                          verbose = 0))
  expect_length(o0[grepl("ignored", o0)], 0)
  o2 <- capture.output(r <- gl.report.maf(testset.gl, plot.display = FALSE,
                                          verbose = 2))
  drop.lines <- o2[grepl("ignored", o2)]
  expect_gt(length(drop.lines), 0)
  expect_true(any(grepl("EmmacNormLeic", drop.lines)))
})
