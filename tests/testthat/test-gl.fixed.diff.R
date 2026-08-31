# Characterization tests for gl.fixed.diff
# Baseline snapshotted before review (review-gl.fixed.diff).
# Assertions marked [approved diff] were flipped in Phase C.

gl4 <- gl.keep.pop(testset.gl,
                   pop.list = c("EmmacBrisWive", "EmmacBurdMist",
                                "EmmacClarJack", "EmmacRussEube"),
                   verbose = 0)

test_that("fd counts, silence, visibility, input untouched", {
  xcopy <- gl4
  o <- capture.output(v <- withVisible(gl.fixed.diff(gl4, verbose = 0)))
  expect_length(o, 0)
  expect_true(v$visible)
  expect_s3_class(v$value, "fd")
  expect_identical(xcopy, gl4)
  fd <- as.matrix(v$value$fd)
  # verified against an independent recomputation from the genotype matrix
  expect_equal(fd["EmmacBrisWive", "EmmacBurdMist"], 1)
  expect_equal(fd["EmmacBrisWive", "EmmacClarJack"], 0)
  expect_equal(fd["EmmacBrisWive", "EmmacRussEube"], 4)
  expect_equal(fd["EmmacBurdMist", "EmmacClarJack"], 1)
  expect_equal(fd["EmmacBurdMist", "EmmacRussEube"], 1)
  expect_equal(fd["EmmacClarJack", "EmmacRussEube"], 2)
  expect_true(isSymmetric(fd))
})

test_that("monomorph handling and the mono.rm parameter", {
  # [approved diff F1] baseline: mono.rm was documented ([default TRUE]
  # removes monomorphs) but never referenced in the body, and the flag
  # logic was inverted - gl.filter.monomorphs ran only when the flag said
  # the data were already monomorph-free. gl4 has 210 monomorphic loci;
  # all 255 were retained and mono.rm had no effect.
  invisible(capture.output(fd.T <- gl.fixed.diff(gl4, verbose = 0)))
  invisible(capture.output(fd.F <- gl.fixed.diff(gl4, mono.rm = FALSE,
                                                 verbose = 0)))
  expect_equal(nLoc(fd.T$gl), 45)                    # monomorphs removed
  expect_lte(max(fd.T$nloc, na.rm = TRUE), 45)
  expect_equal(nLoc(fd.F$gl), 255)                   # FALSE retains them
  expect_equal(max(fd.F$nloc, na.rm = TRUE), 245)    # the old denominators
  # raw fd counts are blind to monomorphic loci - identical either way
  expect_identical(as.matrix(fd.T$fd), as.matrix(fd.F$fd))
})

test_that("tloc>0 with test=TRUE coerces tloc and warns", {
  # [approved diff F2] baseline: the warning printed at verbose 0.
  o <- capture.output(fd <- gl.fixed.diff(gl4, tloc = 0.05, test = TRUE,
                                          reps = 2, verbose = 0))
  expect_length(o, 0)
  # coercion to tloc=0 means fd counts equal the absolute run
  invisible(capture.output(fd0 <- gl.fixed.diff(gl4, verbose = 0)))
  expect_identical(as.matrix(fd$fd), as.matrix(fd0$fd))
  expect_true(isSymmetric(as.matrix(fd$pval)))
})

test_that("fd-class input is accepted and reproduces the matrix", {
  invisible(capture.output(fd1 <- gl.fixed.diff(gl4, verbose = 0)))
  invisible(capture.output(fd2 <- gl.fixed.diff(fd1, verbose = 0)))
  expect_s3_class(fd2, "fd")
  expect_identical(as.matrix(fd2$fd), as.matrix(fd1$fd))
})

test_that("tloc out of range is fatal; fewer than two populations fatal", {
  expect_error(gl.fixed.diff(gl4, tloc = 0.7, verbose = 0))
  expect_error(gl.fixed.diff(gl4, tloc = -0.1, verbose = 0))
  one <- gl.keep.pop(testset.gl, pop.list = "EmmacRussEube", verbose = 0)
  expect_error(suppressWarnings(gl.fixed.diff(one, verbose = 0)))
})

test_that("identical split populations give zero fixed differences", {
  gsub <- gl.keep.pop(testset.gl, pop.list = "EmmacRussEube", verbose = 0)
  pop(gsub) <- factor(rep(c("A", "B"), length.out = nInd(gsub)))
  invisible(capture.output(fd <- gl.fixed.diff(gsub, verbose = 0)))
  expect_equal(as.matrix(fd$fd)["A", "B"], 0)
})
