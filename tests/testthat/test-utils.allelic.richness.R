# Characterization tests for utils.allelic.richness
# Baseline snapshotted before review (kernel-wave review, dev at ddaed27).
# The rarefaction arithmetic was validated against an independent
# recomputation in the gl.report.allelerich review (PR #286).

test_that("the rarefied richness kernel is stable on a fixture", {
  sub <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[c(1, 2, 11)],
                     verbose = 0)
  sub <- gl.filter.allna(sub, verbose = 0)
  gm <- as.matrix(sub[pop(sub) == popNames(sub)[1], ])
  df <- t(gm)
  r <- utils.allelic.richness(df, seq_len(ncol(df)), boot_method = "loc")
  expect_equal(r, 1.2, tolerance = 1e-6)
  # monomorphic-only input: richness 1 per site
  mono <- matrix(0, nrow = 4, ncol = 5)
  r2 <- utils.allelic.richness(t(mono), 1:4, boot_method = "loc")
  expect_equal(r2, 1)
})
