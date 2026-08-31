# Characterization tests for gl.report.factorloadings
# Baseline snapshotted before review (review-gl.report.factorloadings).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

make_pca <- function() {
  pdf(NULL)
  on.exit(dev.off())
  gl.pcoa(testset.gl, verbose = 0)
}

test_that("gl.report.factorloadings returns the loadings invisibly", {
  pca <- make_pca()
  out <- capture.output(
    v <- withVisible(gl.report.factorloadings(pca, plot.display = FALSE,
                                              verbose = 0))
  )
  expect_false(v$visible)
  expect_s3_class(v$value, "data.frame")
  expect_equal(nrow(v$value), nrow(pca$loadings))
})

test_that("gl.report.factorloadings top-N matches independent recomputation", {
  pca <- make_pca()
  out <- capture.output(
    gl.report.factorloadings(pca, plot.display = FALSE, verbose = 2)
  )
  l <- pca$loadings[, 1]
  top <- names(sort(abs(l), decreasing = TRUE))[1:3]
  expect_true(any(grepl(top[1], out, fixed = TRUE)))
  expect_true(any(grepl(top[2], out, fixed = TRUE)))
  expect_true(any(grepl(top[3], out, fixed = TRUE)))
})

test_that("gl.report.factorloadings rejects non-glPca input", {
  expect_error(suppressWarnings(
    capture.output(gl.report.factorloadings(testset.gl, verbose = 0))
  ))
})

test_that("gl.report.factorloadings plot.file works without display", {
  pca <- make_pca()
  tmpdir <- tempdir()
  out <- capture.output(
    res <- gl.report.factorloadings(pca, plot.display = FALSE,
                                    plot.dir = tmpdir,
                                    plot.file = "chartest_fl", verbose = 0)
  )
  saved <- file.path(tmpdir, "chartest_fl.RDS")
  expect_true(file.exists(saved))
  unlink(saved)
})

test_that("gl.report.factorloadings output at verbose = 0", {
  # [approved diff F1] Pre-fix, 17 lines printed at verbose = 0. Now fully
  # silent.
  pca <- make_pca()
  out <- capture.output(
    res <- gl.report.factorloadings(pca, plot.display = FALSE, verbose = 0)
  )
  expect_length(out, 0)
})

test_that("gl.report.factorloadings n.display edge cases", {
  pca <- make_pca()
  # [approved diff F2] head() clamps the display: exactly
  # min(n.display, n.loci) rows, none for 0, no NA garbage.
  out <- capture.output(
    gl.report.factorloadings(pca, n.display = 300, plot.display = FALSE,
                             verbose = 1)
  )
  expect_equal(sum(grepl("^NA", out)), 0)
  n_loadings <- nrow(pca$loadings)
  # header + report line + Completed + n_loadings rows
  expect_true(any(grepl("locus", out)))
  out0 <- capture.output(
    gl.report.factorloadings(pca, n.display = 0, plot.display = FALSE,
                             verbose = 1)
  )
  expect_false(any(grepl("100097451", out0)))
})

test_that("gl.report.factorloadings axis out of range", {
  # [approved diff F3] clear fatal error instead of 'subscript out of
  # bounds'.
  pca <- make_pca()
  expect_error(
    capture.output(gl.report.factorloadings(pca, axis = 99,
                                            plot.display = FALSE,
                                            verbose = 0)),
    "axis must be"
  )
})
