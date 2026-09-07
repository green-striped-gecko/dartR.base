# Characterization tests for gl.save
# Baseline snapshotted before review (review-gl.save).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("saves a loadable file; invisible return", {
  f <- file.path(tempdir(), "test-gl-save.rds")
  o <- capture.output(v <- withVisible(g <- gl.save(testset.gl, f,
                                                    verbose = 0)))
  expect_false(v$visible)
  expect_true(file.exists(f))
  x <- readRDS(f)
  expect_identical(as.matrix(x), as.matrix(testset.gl))
  expect_equal(nInd(x), nInd(testset.gl))
  # [approved diff F1] baseline: two messages printed at verbose 0.
  expect_length(o, 0)
  unlink(f)
})

test_that("returned object is the input", {
  f <- file.path(tempdir(), "test-gl-save2.rds")
  invisible(capture.output(g <- gl.save(testset.gl, f, verbose = 0)))
  # [approved diff F2] baseline: the class-attribute stamping modified
  # the RETURNED object too; now only the saved copy is stamped.
  expect_identical(g, testset.gl)
  unlink(f)
})

test_that("bad path errors", {
  # [approved diff F3] baseline: raw connection error.
  expect_error(suppressWarnings(capture.output(
    gl.save(testset.gl, "Z:/no/such/dir/x.rds", verbose = 0))),
    "does not exist")
})
