# Characterization tests for gl.load
# Baseline snapshotted before review (review-gl.load).
# Assertions marked [approved diff] were flipped in Phase C.

f <- file.path(tempdir(), "test-gl-load.rds")
invisible(capture.output(gl.save(testset.gl, f, verbose = 0)))

test_that("round trip preserves the data; invisible return", {
  o <- capture.output(v <- withVisible(x1 <- gl.load(f, verbose = 0)))
  expect_false(v$visible)
  expect_identical(as.matrix(v$value), as.matrix(testset.gl))
  expect_equal(nInd(v$value), nInd(testset.gl))
  expect_equal(nLoc(v$value), nLoc(testset.gl))
  expect_equal(nPop(v$value), nPop(testset.gl))
})

test_that("verbose 0 silence and the compliance parameter", {
  # BASELINE: 21 lines print at verbose 0 (the ungated 'Loaded object'
  # message plus the full gl.compliance.check chatter), and
  # compliance = TRUE vs FALSE are byte-identical - the parameter is
  # inert because the check runs unconditionally.
  # [approved diff F1+F2] baseline: 21 lines printed at verbose 0 and
  # compliance TRUE/FALSE were byte-identical (the check ran
  # unconditionally).
  o <- capture.output(x <- gl.load(f, verbose = 0))
  expect_length(o, 0)
  ot <- capture.output(xt <- gl.load(f, compliance = TRUE, verbose = 0))
  of <- capture.output(xf <- gl.load(f, compliance = FALSE, verbose = 0))
  expect_length(ot, 0)
  expect_length(of, 0)
})

test_that("compliance visible only when requested at verbose 2", {
  # [approved diff F1] baseline: the compliance banner showed for BOTH
  # settings.
  ot <- capture.output(xt <- gl.load(f, compliance = TRUE, verbose = 2))
  of <- capture.output(xf <- gl.load(f, compliance = FALSE, verbose = 2))
  expect_true(any(grepl("compliance", ot)))
  expect_false(any(grepl("compliance", of)))
})

test_that("missing and non-genlight files", {
  # [approved diff F3] baseline: raw connection error / cryptic
  # "no applicable method for `@`" error.
  expect_error(suppressWarnings(capture.output(
    gl.load("no-such-file-xyz.rds", verbose = 0))), "not found")
  f2 <- file.path(tempdir(), "test-notgl.rds")
  saveRDS(data.frame(a = 1:3), f2)
  expect_error(capture.output(gl.load(f2, verbose = 0)),
               "genlight")
  unlink(f2)
})

test_that("dartR objects saved by another package (legacy dartR) load without loading it", {
  # An object of class "dartR" defined by another package, without the fbm
  # slot, as saved by the legacy dartR; the package it points to is not
  # installed, so any class lookup through it would fail
  x <- methods::as(testset.gl[1:10, 1:20], "genlight")
  attr(x, "class") <- structure("dartR", package = "dartRlegacyNotInstalled")
  expect_false(methods::.hasSlot(x, "fbm"))
  f3 <- file.path(tempdir(), "test-legacy-dartR.rds")
  saveRDS(x, f3)
  on.exit(unlink(f3), add = TRUE)

  o <- capture.output(y <- gl.load(f3, verbose = 0))
  expect_length(o, 0)
  expect_false("dartRlegacyNotInstalled" %in% loadedNamespaces())
  expect_identical(attr(class(y), "package"), "dartR.base")
  expect_true(methods::.hasSlot(y, "fbm"))
  expect_true(methods::validObject(y))
  expect_identical(unname(as.matrix(y)),
                   unname(as.matrix(testset.gl[1:10, 1:20])))
  expect_identical(y@other$ind.metrics, testset.gl[1:10, 1:20]@other$ind.metrics)
  expect_true(any(grepl("another package",
                        capture.output(gl.load(f3, verbose = 2)))))
  # and to FBM on request, silently
  o <- capture.output(z <- gl.load(f3, fbm = TRUE, verbose = 0))
  expect_length(o, 0)
  expect_false(is.null(z@fbm))
})

test_that("an FBM object saved with saveRDS whose backing file is gone errors clearly", {
  skip_if_not_installed("bigstatsr")
  x <- gl.gen2fbm(testset.gl[1:5, 1:10], verbose = 0)
  f4 <- file.path(tempdir(), "test-raw-fbm.rds")
  saveRDS(x, f4)
  on.exit(unlink(f4), add = TRUE)
  # while the backing file exists it loads
  expect_identical(as.matrix(gl.load(f4, verbose = 0)), as.matrix(x))
  bk <- x@fbm$backingfile
  rm(x); invisible(gc())
  unlink(bk)
  expect_error(capture.output(gl.load(f4, verbose = 0)), "gl.save")
})
