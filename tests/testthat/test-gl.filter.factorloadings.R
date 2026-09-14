# Characterization tests for gl.filter.factorloadings
# Baseline snapshotted before review (review-gl.filter.factorloadings).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

prep_pca <- function() {
  gl.pcoa(testset.gl, verbose = 0)
}

test_that("gl.filter.factorloadings retain semantics", {
  pca <- prep_pca()
  expect_equal(nrow(pca$loadings), 111)  # loci after monomorph removal
  out1 <- capture.output(
    f_false <- gl.filter.factorloadings(testset.gl, pca = pca,
                                        threshold = 0.1, retain = FALSE,
                                        plot.display = FALSE, verbose = 0)
  )
  expect_equal(nLoc(f_false), 80)
  expect_equal(nrow(f_false@other$loc.metrics), nLoc(f_false))
  expect_true(all(abs(pca$loadings[match(locNames(f_false),
                                         locNames(gl.filter.monomorphs(
                                           testset.gl, verbose = 0))), 1])
                  < 0.1))
  # [approved diff F1] retain = TRUE now keeps the high-loading loci as
  # documented (previously it returned the identical 80-locus complement)
  f_true <- gl.filter.factorloadings(testset.gl, pca = pca,
                                     threshold = 0.1, retain = TRUE,
                                     plot.display = FALSE, verbose = 0)
  expect_equal(nLoc(f_true), 31)
  xm <- gl.filter.monomorphs(testset.gl, verbose = 0)
  hi <- locNames(xm)[abs(pca$loadings[, 1]) >= 0.1]
  expect_setequal(locNames(f_true), hi)
  # the two retain values partition the polymorphic loci
  expect_setequal(c(locNames(f_true), locNames(f_false)), locNames(xm))
})

test_that("gl.filter.factorloadings rejects a mismatched pca", {
  # [approved diff F2] a pca from a different locus set previously produced
  # a recycling warning and misaligned filtering; it now errors clearly
  pca <- prep_pca()
  xs <- testset.gl[, 1:150]
  xs@other$loc.metrics <- testset.gl@other$loc.metrics[1:150, ]
  expect_error(
    gl.filter.factorloadings(xs, pca = pca, threshold = 0.1,
                             plot.display = FALSE, verbose = 0),
    "polymorphic loci"
  )
})

test_that("gl.filter.factorloadings appends a single history entry", {
  # [approved diff F3] previously two internal entries
  # (gl.filter.monomorphs, gl.keep.loc) and none for the function itself
  pca <- prep_pca()
  h0 <- length(testset.gl@other$history)
  f <- gl.filter.factorloadings(testset.gl, pca = pca, threshold = 0.1,
                                plot.display = FALSE, verbose = 0)
  expect_equal(length(f@other$history), h0 + 1)
  expect_true(grepl("gl.filter.factorloadings",
                    deparse(f@other$history[[h0 + 1]])[1]))
})

test_that("gl.filter.factorloadings is silent at verbose 0 and input untouched", {
  pca <- prep_pca()
  out <- capture.output(
    f <- gl.filter.factorloadings(testset.gl, pca = pca, threshold = 0.1,
                                  verbose = 0)
  )
  expect_length(out, 0)
  expect_identical(testset.gl, dartR.data::testset.gl)
})

test_that("gl.filter.factorloadings saves the plot with plot.display = FALSE", {
  pca <- prep_pca()
  pd <- file.path(tempdir(), "fl-plots")
  dir.create(pd, showWarnings = FALSE)
  expect_no_error(
    gl.filter.factorloadings(testset.gl, pca = pca, threshold = 0.1,
                             plot.display = FALSE, plot.file = "fltest",
                             plot.dir = pd, verbose = 0)
  )
  expect_true(length(list.files(pd, pattern = "fltest")) > 0)
})

test_that("gl.filter.factorloadings error paths", {
  pca <- prep_pca()
  # [approved diff F5] axis out of range now names the valid range
  expect_error(
    gl.filter.factorloadings(testset.gl, pca = pca, threshold = 0.1,
                             axis = 99, plot.display = FALSE, verbose = 0),
    "axis must be between"
  )
  # [approved diff F4] non-glPca objects raise a labelled error
  expect_error(
    gl.filter.factorloadings(testset.gl, pca = testset.gl, threshold = 0.1,
                             plot.display = FALSE, verbose = 0),
    "glPca"
  )
  # [approved diff F7] missing threshold raises a labelled error
  expect_error(
    gl.filter.factorloadings(testset.gl, pca = pca,
                             plot.display = FALSE, verbose = 0),
    "threshold is required"
  )
})
