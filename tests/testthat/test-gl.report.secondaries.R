# Characterization tests for gl.report.secondaries
# Baseline snapshotted before review (review-gl.report.secondaries).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.report.secondaries returns the documented 8-parameter data.frame", {
  skip_if_not_installed("dartR.data")
  library(dartR.data)
  res <- suppressWarnings(gl.report.secondaries(platypus.gl, verbose = 0))
  expect_s3_class(res, "data.frame")
  expect_equal(dim(res), c(8, 2))
  expect_equal(
    res$Param,
    c("n.total.tags", "n.SNPs.secondaries", "n.invariant.tags",
      "n.tags.secondaries", "n.inv.gen", "mean.len.tag", "n.invariant",
      "Lambda")
  )
})

test_that("gl.report.secondaries values match independent recomputation (platypus.gl)", {
  skip_if_not_installed("dartR.data")
  library(dartR.data)
  res <- suppressWarnings(gl.report.secondaries(platypus.gl, verbose = 0))
  vals <- res$Value
  names(vals) <- res$Param
  # Independently recomputed from AlleleID / TrimmedSequence (see review report)
  expect_equal(unname(vals["n.total.tags"]), 991)
  expect_equal(unname(vals["n.SNPs.secondaries"]), 9)
  expect_equal(unname(vals["n.invariant.tags"]), 51288)
  expect_equal(unname(vals["n.tags.secondaries"]), 9)
  expect_equal(unname(vals["n.inv.gen"]), 65820)
  expect_equal(unname(vals["mean.len.tag"]), 67.42684, tolerance = 1e-6)
  expect_equal(unname(vals["n.invariant"]), 3524008)
  expect_equal(unname(vals["Lambda"]), 0.01913782, tolerance = 1e-6)
})

test_that("gl.report.secondaries leaves the input object untouched (report family)", {
  skip_if_not_installed("dartR.data")
  library(dartR.data)
  x <- platypus.gl
  h_before <- length(x@other$history)
  invisible(suppressWarnings(gl.report.secondaries(x, verbose = 0)))
  expect_identical(x, platypus.gl)
  expect_equal(length(x@other$history), h_before)
})

test_that("gl.report.secondaries rejects SilicoDArT data", {
  expect_error(gl.report.secondaries(testset.gs, verbose = 0))
})

test_that("gl.report.secondaries on a dataset with no secondaries", {
  # [approved diff F1] Pre-fix this crashed ('Subsetting resulted in zero
  # loci') via an unused x[, duplicated(b)] subset. Post-fix the documented
  # no-secondaries branch runs: returns the 8-row data.frame with
  # secondaries-related values 0/NA. testset.gl has 255 loci, all on
  # unique tags.
  res <- gl.report.secondaries(testset.gl, verbose = 0)
  expect_s3_class(res, "data.frame")
  expect_equal(dim(res), c(8, 2))
  vals <- res$Value
  names(vals) <- res$Param
  expect_equal(unname(vals["n.total.tags"]), 255)
  expect_equal(unname(vals["n.SNPs.secondaries"]), 0)
  expect_true(is.na(vals["n.invariant.tags"]))
  expect_true(is.na(vals["Lambda"]))
})

test_that("gl.report.secondaries verbosity gating at verbose = 0", {
  skip_if_not_installed("dartR.data")
  library(dartR.data)
  # [approved diff F2/F3] Pre-fix the results block (11 lines) and the
  # TrimmedSequence warning printed unconditionally; both now gate at
  # verbose >= 1, so verbose = 0 is fully silent. The result is assigned
  # inside capture.output because the visible return (kept by decision,
  # review F5 rejected) would otherwise auto-print.
  out <- capture.output(suppressWarnings(
    res <- gl.report.secondaries(platypus.gl, verbose = 0)
  ))
  expect_length(out, 0)
  xx <- platypus.gl
  xx@other$loc.metrics$TrimmedSequence <- NULL
  out2 <- capture.output(suppressWarnings(
    res2 <- gl.report.secondaries(xx, verbose = 0)
  ))
  expect_length(out2, 0)
})

test_that("gl.report.secondaries returns visibly", {
  skip_if_not_installed("dartR.data")
  library(dartR.data)
  # Visible return retained by decision (review F5 rejected).
  v <- withVisible(suppressWarnings(
    gl.report.secondaries(platypus.gl, verbose = 0)
  ))
  expect_true(v$visible)
})
