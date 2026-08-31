# Characterization tests for gl.Ho and gl.He
# Baseline snapshotted before review (review-gl.Ho-gl.He).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("gl.Ho matches hand computation, named, full length", {
  m <- as.matrix(testset.gl)
  ho <- gl.Ho(testset.gl)
  expect_identical(ho, colMeans(m == 1, na.rm = TRUE))
  expect_length(ho, nLoc(testset.gl))
  expect_identical(names(ho), locNames(testset.gl))
})

test_that("gl.He matches hand computation and gl.alf", {
  m <- as.matrix(testset.gl)
  he <- gl.He(testset.gl)
  p <- colMeans(m, na.rm = TRUE) / 2
  expect_identical(he, 2 * p * (1 - p))
  alf <- gl.alf(testset.gl)
  expect_equal(unname(he), 2 * alf$alf1 * alf$alf2)
})

test_that("per-locus Ho aggregates to the report's population Ho", {
  pdf(NULL); on.exit(dev.off())
  sub <- testset.gl[as.character(pop(testset.gl)) == "EmmacBurdMist", ]
  r <- as.data.frame(gl.report.heterozygosity(testset.gl,
                                              plot.display = FALSE,
                                              verbose = 0))
  expect_equal(round(mean(gl.Ho(sub), na.rm = TRUE), 6),
               r[r$pop == "EmmacBurdMist", "Ho"])
})

test_that("pure-function behaviour: silent, visible vector, input untouched", {
  o <- capture.output(v <- withVisible(gl.Ho(testset.gl)))
  expect_length(o, 0)
  expect_true(v$visible)  # the vector IS the product; visible is correct
  xcopy <- testset.gl
  invisible(gl.Ho(xcopy)); invisible(gl.He(xcopy))
  expect_identical(xcopy, testset.gl)
})

test_that("all-NA loci yield NaN", {
  m2 <- matrix(c(0,1,2, NA,NA,NA, 2,1,0), ncol = 3,
               dimnames = list(paste0("I",1:3), c("L1","LNA","L3")))
  gx <- new("genlight", gen = m2, ploidy = 2)
  expect_true(is.nan(gl.Ho(gx)["LNA"]))
  expect_true(is.nan(gl.He(gx)["LNA"]))
})

test_that("SilicoDArT input", {
  # [approved diff F1] presence/absence data is now a fatal error instead
  # of silently returning meaningless values.
  expect_error(gl.Ho(testset.gs))
  expect_error(gl.He(testset.gs))
})
