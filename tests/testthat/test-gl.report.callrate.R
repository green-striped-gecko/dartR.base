test_that("gl.report.callrate returns the input identical and unaltered (SNP, both methods)", {
  x <- testset.gl
  out <- capture.output(x2 <- gl.report.callrate(x, plot.display = FALSE, verbose = 0))
  expect_identical(x2, x)
  out2 <- capture.output(x3 <- gl.report.callrate(x, method = "ind", plot.display = FALSE, verbose = 0))
  expect_identical(x3, x)
})

test_that("gl.report.callrate works on SilicoDArT", {
  x <- testset.gs
  out <- capture.output(x2 <- gl.report.callrate(x, plot.display = FALSE, verbose = 0))
  expect_true(inherits(x2, "genlight"))
})

test_that("reported quantile table is consistent with an independent computation", {
  x <- testset.gl
  x <- utils.recalc.callrate(x, verbose = 0)
  cr <- x@other$loc.metrics$CallRate
  out <- capture.output(gl.report.callrate(x, plot.display = FALSE, verbose = 2))
  med_line <- grep("Median", out, value = TRUE)
  expect_length(med_line, 1)
  reported_med <- as.numeric(sub(".*:\\s*", "", med_line))
  expect_equal(reported_med, as.numeric(summary(cr)[3]), tolerance = 1e-4)
})

test_that("verbose = 0 is fully silent; verbose >= 1 prints the report (F1 fix, VRB5)", {
  x <- testset.gl
  out0 <- capture.output(gl.report.callrate(x, plot.display = FALSE, verbose = 0))
  expect_length(out0, 0)
  out0i <- capture.output(gl.report.callrate(x, method = "ind", plot.display = FALSE, verbose = 0))
  expect_length(out0i, 0)
  out1 <- capture.output(gl.report.callrate(x, plot.display = FALSE, verbose = 1))
  expect_true(any(grepl("Reporting Call Rate by Locus", out1)))
})

test_that("unknown method warns and coerces to 'loc' (F3 fix)", {
  x <- testset.gl
  out <- capture.output(rp <- gl.report.callrate(x, method = "pop", plot.display = FALSE, verbose = 1))
  expect_true(any(grepl("method must be either", out)))
  expect_true(any(grepl("Reporting Call Rate by Locus", out)))
  expect_identical(rp, x)
})

test_that("stale CallRate comes back unaltered (F4 fix); report still uses fresh values", {
  x <- testset.gl
  x_stale <- x
  x_stale@other$loc.metrics$CallRate <- rep(0.5, nLoc(x))
  out <- capture.output(x_back <- gl.report.callrate(x_stale, plot.display = FALSE, verbose = 1))
  expect_identical(x_back, x_stale)
  # the report used recalculated values, not the planted 0.5s
  med_line <- grep("Median", out, value = TRUE)
  expect_false(grepl(": *0\\.5 *$", med_line))
})

test_that("no stray ')' in output; ind.to.list = 0 skips the individual listing (F5/F7 fixes)", {
  x <- testset.gl
  outi <- capture.output(gl.report.callrate(x, method = "ind", plot.display = FALSE, verbose = 2))
  expect_false(any(grepl("^\\)", outi)))
  out_zero <- capture.output(gl.report.callrate(x, method = "ind", ind.to.list = 0, plot.display = FALSE, verbose = 2))
  expect_false(any(grepl("individuals with the lowest CallRates", out_zero)))
})
