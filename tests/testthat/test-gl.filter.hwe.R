# Characterization tests for gl.filter.hwe (post excess.het migration)
# Baseline snapshotted on branch migrate-excess.het-to-hwe. Assertions
# marked [approved diff] were flipped in Phase C to reflect approved
# behaviour changes.

test_that("gl.filter.hwe direction='excess' removes exactly the report's loci", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  e <- suppressWarnings(gl.report.hwe(testset.gl, method_sig = "ChiSquare",
                                      multi_comp = TRUE,
                                      multi_comp_method = "fdr",
                                      direction = "excess",
                                      plot.out = FALSE, verbose = 0))
  f <- suppressWarnings(gl.filter.hwe(testset.gl, test.type = "ChiSquare",
                                      mult.comp.adj = TRUE,
                                      mult.comp.adj.method = "fdr",
                                      direction = "excess", verbose = 0))
  removed <- setdiff(locNames(testset.gl), locNames(f))
  expect_setequal(removed, unique(as.character(e$Locus)))
  expect_equal(nrow(f@other$loc.metrics), nLoc(f))
})

test_that("gl.filter.hwe default direction='both' is the pre-migration behaviour", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  f <- suppressWarnings(gl.filter.hwe(testset.gl, verbose = 0))
  expect_equal(nLoc(f), 249)
  expect_equal(nrow(f@other$loc.metrics), nLoc(f))
})

test_that("gl.filter.hwe appends one history entry, invisible return, input untouched", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  h <- length(testset.gl@other$history)
  v <- suppressWarnings(withVisible(gl.filter.hwe(testset.gl, verbose = 0)))
  expect_false(v$visible)
  expect_equal(length(v$value@other$history), h + 1)
  xcopy <- testset.gl
  suppressWarnings(invisible(gl.filter.hwe(xcopy, verbose = 0)))
  expect_identical(xcopy, testset.gl)
})

test_that("gl.filter.excess.het stub matches the original's removals on LBP", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  expect_warning(
    f <- gl.filter.excess.het(LBP, Yates = TRUE, verbose = 0),
    "deprecated"
  )
  # the original implementation removed 6 loci from LBP with Yates=TRUE
  expect_equal(nLoc(LBP) - nLoc(f), 6)
  expect_equal(length(f@other$history) - length(LBP@other$history), 1)
})

test_that("gl.filter.hwe is silent at verbose = 0", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  out <- capture.output(
    f <- suppressWarnings(gl.filter.hwe(testset.gl, verbose = 0))
  )
  expect_length(out, 0)
})
