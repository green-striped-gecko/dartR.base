# Characterization tests for gl.report.hwe (post excess.het migration)
# Baseline snapshotted on branch migrate-excess.het-to-hwe. Assertions
# marked [approved diff] were flipped in Phase C to reflect approved
# behaviour changes.

test_that("gl.report.hwe direction partitions the significant set", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  b <- suppressWarnings(gl.report.hwe(testset.gl, method_sig = "ChiSquare",
                                      direction = "both",
                                      plot.out = FALSE, verbose = 0))
  e <- suppressWarnings(gl.report.hwe(testset.gl, method_sig = "ChiSquare",
                                      direction = "excess",
                                      plot.out = FALSE, verbose = 0))
  d <- suppressWarnings(gl.report.hwe(testset.gl, method_sig = "ChiSquare",
                                      direction = "deficit",
                                      plot.out = FALSE, verbose = 0))
  expect_true(all(e$Het > e$Het.exp))
  expect_true(all(d$Het < d$Het.exp))
  expect_equal(nrow(e) + nrow(d) + sum(b$Het == b$Het.exp), nrow(b))
})

test_that("gl.report.hwe Het.exp matches independent recomputation", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  df <- suppressWarnings(gl.report.hwe(testset.gl, sig_only = FALSE,
                                       plot.out = FALSE, verbose = 0))
  df <- as.data.frame(df)
  p <- (2 * df$Hom_1 + df$Het) / (2 * df$N)
  expect_equal(df$Het.exp, 2 * df$N * p * (1 - p))
})

test_that("gl.report.hwe min.hobs screens the tested set", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  df <- suppressWarnings(gl.report.hwe(testset.gl, sig_only = FALSE,
                                       min.hobs = 0.5,
                                       plot.out = FALSE, verbose = 0))
  df <- as.data.frame(df)
  expect_true(all(df$Het / df$N >= 0.5))
})

test_that("gl.report.hwe returns invisibly and leaves the input untouched", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  v <- suppressWarnings(withVisible(gl.report.hwe(testset.gl,
                                                  plot.out = FALSE,
                                                  verbose = 0)))
  expect_false(v$visible)
  xcopy <- testset.gl
  suppressWarnings(invisible(gl.report.hwe(xcopy, plot.out = FALSE,
                                           verbose = 0)))
  expect_identical(xcopy, testset.gl)
})

test_that("gl.report.hwe chi-square p-values match HWChisq directly", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  df <- suppressWarnings(gl.report.hwe(testset.gl, method_sig = "ChiSquare",
                                       sig_only = FALSE,
                                       plot.out = FALSE, verbose = 0))
  r1 <- as.data.frame(df)[1, ]
  pv <- suppressWarnings(HardyWeinberg::HWChisq(
    c(AA = r1$Hom_1, AB = r1$Het, BB = r1$Hom_2), verbose = FALSE)$pval)
  expect_equal(r1$Prob, pv)
})

test_that("gl.report.hwe tested set no longer depends on verbosity", {
  # [approved diff R1] Pre-fix, the skipping of monomorphic and small
  # populations sat inside verbose >= 2 gates (verbose 0: 311 rows across
  # 30 pops; verbose 2: 245 rows across 23). Post-fix both verbosities
  # test the same 23 qualifying populations.
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  a0 <- suppressWarnings(gl.report.hwe(testset.gl, sig_only = FALSE,
                                       plot.out = FALSE, verbose = 0))
  invisible(capture.output(
    a2 <- suppressWarnings(gl.report.hwe(testset.gl, sig_only = FALSE,
                                         plot.out = FALSE, verbose = 2))
  ))
  # 196 rows with dartR.data 1.2.2 (all-monomorphic pops now skipped
  # instead of crashing in gl.filter.monomorphs)
  expect_equal(nrow(a0), 196)
  expect_identical(as.data.frame(a0), as.data.frame(a2))
})

test_that("gl.report.excess.het stub reproduces the published workflow shape", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  expect_warning(
    res <- gl.report.excess.het(LBP, Yates = TRUE, plot.display = FALSE,
                                verbose = 0),
    "deprecated"
  )
  expect_named(res, c("results.table", "removed.loci"))
  # LBP fixture: the original implementation flagged these 6 loci
  expect_equal(length(res$removed.loci), 6)
})
