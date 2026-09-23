# Characterization baseline (function-review, 2026-09-24): outputs of the
# reviewed state. Detects change; does not assert correctness.
test_that("baseline: hwe_summary on bandicoot.gl", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  capture.output(r <- suppressWarnings(
    gl.diagnostics.hwe(bandicoot.gl, stdErr = FALSE, verbose = 0)))
  s <- as.data.frame(r$hwe_summary)
  expect_equal(s$Population, c("NSW", "QLD", "SA", "VIC", "WA"))
  expect_equal(s$nSig, c(22, 27, 53, 25, 20))
  # [approved diff, change 4] alpha x tests run per population, not x nLoc
  expect_equal(s$nExpected, 0.05 * c(976, 953, 971, 955, 969))
  expect_equal(s$Deficiency, c(18, 24, 52, 24, 14))
  expect_equal(s$Excess, c(4, 3, 1, 1, 6))
  expect_equal(s$ChiSquare, c(1151.773124, 1158.936715, 1541.817417,
                              1055.195549, 1252.109628), tolerance = 1e-8)
  expect_equal(s$pvalue, rep(1, 5))
  expect_named(r, "hwe_summary")
})

test_that("jackknife standard errors on 100 loci", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  capture.output(r <- suppressWarnings(
    gl.diagnostics.hwe(bandicoot.gl[, 1:100], stdErr = TRUE, n.cores = 1,
                       verbose = 0)))
  expect_named(r, c("hwe_summary", "StdErr"))
  # [approved diff, changes 1-2] reviewed state: 1.408998498e-04,
  # 2.937462658e-05 (rounded leave-one-out values, SE without the n - 1
  # factor). Recompute by brute force with the standard jackknife SE.
  x <- bandicoot.gl[, 1:100]
  L <- nLoc(x)
  loo <- t(vapply(seq_len(L), function(i) {
    o <- utils.basic.stats(gl.drop.loc(x, locNames(x)[i], verbose = 0),
                           rounded = FALSE)$overall
    c(o[["Fis"]], o[["Fst"]])
  }, numeric(2)))
  se <- apply(loo, 2, function(v) sqrt((L - 1) / L * sum((v - mean(v))^2)))
  expect_equal(unname(r$StdErr), se, tolerance = 1e-10)
})

test_that("the barplot counts loci in every bar, against the expected null", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  d <- tempfile(); dir.create(d)
  capture.output(suppressWarnings(
    gl.diagnostics.hwe(bandicoot.gl, stdErr = FALSE, plot.file = "p",
                       plot.dir = d, verbose = 0)))
  p <- readRDS(file.path(d, "p.RDS"))
  bars <- p[[2]]$data
  obs <- bars[bars$Data == "Observed", ]
  null <- bars[bars$Data == "Null expectation", ]
  # 991 loci tested; the observed bars and the null both sum to it
  expect_equal(sum(obs$Freq), 991)
  expect_equal(sum(null$Freq), 991, tolerance = 1e-6)
  expect_equal(obs$Freq[obs$nPop == "0"], 863)
  # histogram line at the uniform expectation
  h <- p[[1]]
  expect_equal(ggplot2::layer_data(h, 2)$yintercept[1], 4824 / 20)
})

test_that("gl.diagnostics.hwe is silent at verbose = 0", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  out <- capture.output(suppressWarnings(
    gl.diagnostics.hwe(bandicoot.gl[, 1:100], verbose = 0)), type = "output")
  msg <- capture.output(suppressWarnings(
    gl.diagnostics.hwe(bandicoot.gl[, 1:100], verbose = 0)), type = "message")
  expect_length(out, 0)
  expect_length(msg, 0)
})
