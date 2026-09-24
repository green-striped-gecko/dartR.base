# Characterization baseline (function-review, 2026-09-24): outputs of the
# reviewed state. Detects change; does not assert correctness.
test_that("baseline: gen -> FBM -> gen round trip preserves the object", {
  skip_if_not_installed("bigstatsr")
  x <- testset.gl
  capture.output(f <- gl.gen2fbm(x, verbose = 0))
  expect_length(f@gen, 0)
  expect_s4_class(f@fbm, "FBM.code256")
  expect_equal(unname(as.matrix(f)) + 0, unname(as.matrix(x)) + 0)
  capture.output(b <- gl.fbm2gen(f, verbose = 0))
  expect_null(b@fbm)
  expect_identical(as.matrix(b), as.matrix(x))
  for (s in c("ind.names", "loc.names", "loc.all", "ploidy", "pop",
              "position", "chromosome", "n.loc", "other")) {
    expect_identical(slot(b, s), slot(x, s), label = s)
  }
  # small blocks give the same result
  capture.output(f7 <- gl.gen2fbm(x, chunk = 7L, verbose = 0))
  expect_equal(unname(as.matrix(f7)) + 0, unname(as.matrix(x)) + 0)
  # already FBM-backed: returned as is
  capture.output(f2 <- gl.gen2fbm(f, verbose = 0))
  expect_identical(f2@fbm$backingfile, f@fbm$backingfile)
})

test_that("baseline: unsupported inputs", {
  skip_if_not_installed("bigstatsr")
  expect_error(capture.output(gl.gen2fbm(testset.gs, verbose = 0)),
               "Only SNP data")
  # [approved diff, change 2] a plain genlight failed in both functions
  g <- new("genlight", as.matrix(testset.gl)[1:20, 1:30], ploidy = 2)
  capture.output(gf <- gl.gen2fbm(g, verbose = 0))
  expect_s4_class(gf, "dartR")
  expect_equal(unname(as.matrix(gf)) + 0, unname(as.matrix(g)) + 0)
  capture.output(gb <- gl.fbm2gen(g, verbose = 0))
  expect_identical(gb, g)
  # gen-backed input to fbm2gen is returned unchanged
  capture.output(b <- gl.fbm2gen(testset.gl, verbose = 0))
  expect_identical(b, testset.gl)
})

test_that("gl.fbm2gen gives the same object for any block size", {
  skip_if_not_installed("bigstatsr")
  x <- testset.gl
  capture.output(f <- gl.gen2fbm(x, verbose = 0))
  for (k in c(1L, 7L, 249L, 250L, 1000L)) {
    capture.output(b <- gl.fbm2gen(f, chunk = k, verbose = 0))
    expect_identical(as.matrix(b), as.matrix(x), label = paste("chunk", k))
    expect_identical(b@ploidy, x@ploidy)
    expect_null(b@fbm)
  }
})

test_that("an object with both slots populated is reported as invalid", {
  skip_if_not_installed("bigstatsr")
  capture.output(f <- gl.gen2fbm(testset.gl, verbose = 0))
  bad <- f
  bad@gen <- testset.gl@gen
  expect_error(capture.output(gl.gen2fbm(bad, verbose = 0)),
               "both @fbm and @gen")
})
