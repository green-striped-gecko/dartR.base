# Characterization tests for the utils.hwe helpers
# Baseline snapshotted before review (kernel-wave review, dev at ddaed27).

test_that("GenerateSamples enumerates all genotype compositions", {
  gs <- GenerateSamples(4)
  expect_equal(nrow(gs), choose(4 + 2, 2))
  expect_true(all(rowSums(gs) == 4))
  expect_equal(colnames(gs), c("AA", "AB", "BB"))
})

test_that("CritSam returns a critical-sample matrix", {
  skip_if_not_installed("HardyWeinberg")
  cs <- CritSam(n = 10, Dpos = TRUE, alphalimit = 0.05,
                pvaluetype = "selome")
  expect_true(is.matrix(cs$Xn))
  expect_true(all(abs(rowSums(cs$Xn) - 1) < 1e-9))
  cc <- CritSam_Chi(n = 10, Dpos = TRUE, alphalimit = 0.05, cc = 0.5)
  expect_true(is.matrix(cc$Xn))
})
