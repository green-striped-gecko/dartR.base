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

test_that("FBM objects are saved as gen-backed copies and round-trip", {
  skip_if_not_installed("bigstatsr")
  xf <- gl.gen2fbm(testset.gl, verbose = 0)
  bk <- xf@fbm$backingfile
  f5 <- file.path(tempdir(), "test-save-fbm.rds")
  on.exit(unlink(f5), add = TRUE)
  o <- capture.output(ret <- gl.save(xf, f5, verbose = 0))
  expect_length(o, 0)
  # input untouched: still FBM-backed, same backing file
  expect_identical(ret, xf)
  expect_identical(ret@fbm$backingfile, bk)
  # the file holds genotypes (not a reference to the backing file)
  s <- readRDS(f5)
  expect_null(s@fbm)
  expect_length(s@gen, nInd(xf))
  # same size as saving the gen-backed object
  f6 <- file.path(tempdir(), "test-save-gen.rds")
  on.exit(unlink(f6), add = TRUE)
  invisible(capture.output(gl.save(testset.gl, f6, verbose = 0)))
  expect_equal(file.size(f5), file.size(f6))
  # round trip in both flavours
  expect_identical(as.matrix(gl.load(f5, verbose = 0)), as.matrix(xf))
  yf <- gl.load(f5, fbm = TRUE, verbose = 0)
  expect_false(is.null(yf@fbm))
  expect_identical(as.matrix(yf), as.matrix(xf))
})
