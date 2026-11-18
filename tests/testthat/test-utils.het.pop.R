# Characterization tests for utils.het.pop
# Baseline snapshotted before review (kernel-wave review, dev at ddaed27).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("expected heterozygosity matches the documented computation", {
  sub <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[c(1, 2, 11)],
                     verbose = 0)
  sub <- gl.filter.allna(sub, verbose = 0)
  hp <- utils.het.pop(sub)
  gm <- as.matrix(sub[pop(sub) == popNames(sub)[1], ])
  p0 <- colMeans(gm == 0, na.rm = TRUE); p2 <- colMeans(gm == 2, na.rm = TRUE)
  ph <- colMeans(gm == 1, na.rm = TRUE)
  p <- (2 * p0 + ph) / 2; q <- (2 * p2 + ph) / 2
  H <- 1 - (p^2 + q^2)
  nl <- colSums(!is.na(gm))
  mean_n <- mean(nl[nl > 0])
  expect_equal(hp[[1]],
               round(mean((2 * mean_n / (2 * mean_n - 1)) * H, na.rm = TRUE), 6))
  # [approved diff K5] baseline: the vector came back unnamed.
  expect_equal(names(hp), popNames(sub))  # [approved diff K5]
})
