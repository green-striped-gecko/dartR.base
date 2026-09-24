# Characterization tests for utils.pa.Chao
# Baseline snapshotted before review (review-utils.pa.Chao).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("no private alleles gives 0 in both directions", {
  m <- matrix(c(0, 1, 2, 0, 1, 2), ncol = 2)
  expect_equal(utils.pa.Chao(m, m), list(0, 0))
})

test_that("(n - 1) / n uses the number of allele copies sampled", {
  # two private singletons (f1 = 2, f2 = 0): f1 (f1 - 1) / 2 = 1.
  # [approved diff change 1] baseline: scaled by (S - 1) / S with S = 2
  # private alleles, giving 0.5. Now n = 12 allele copies (6 individuals
  # in the pooled pair), giving 11 / 12.
  m1 <- matrix(c(0, 0, 0, 0, 0, 1, 0, 1, 0), ncol = 3)
  m2 <- matrix(0, 3, 3)
  expect_equal(utils.pa.Chao(m1, m2), list(11 / 12, 0))
  # a single private singleton: f1 (f1 - 1) / 2 = 0
  expect_equal(utils.pa.Chao(m1[, 1:2], m2[, 1:2])[[1]], 0)
})

test_that("platypus pairs", {
  pops <- adegenet::seppop(platypus.gl)
  r <- utils.pa.Chao(pops$SEVERN_ABOVE, pops$SEVERN_BELOW)
  # [approved diff change 1] baseline: 62.60081, 21.84557 and 4.925676
  # with n = number of private alleles
  expect_lt(abs(r[[1]] - 62.60081) / 62.60081, 0.02)
  expect_lt(abs(r[[2]] - 21.84557) / 21.84557, 0.02)
  r <- utils.pa.Chao(pops$SEVERN_BELOW, pops$TENTERFIELD)
  expect_lt(abs(r[[1]] - 4.925676) / 4.925676, 0.03)
})

test_that("FBM-backed pairs give the same estimates", {
  f <- gl.gen2fbm(platypus.gl, verbose = 0)
  pf <- adegenet::seppop(f)
  pg <- adegenet::seppop(platypus.gl)
  expect_equal(utils.pa.Chao(pf$SEVERN_ABOVE, pf$TENTERFIELD),
               utils.pa.Chao(pg$SEVERN_ABOVE, pg$TENTERFIELD))
})
