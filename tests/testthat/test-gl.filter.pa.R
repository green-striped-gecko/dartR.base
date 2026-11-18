# Characterization tests for gl.filter.pa
# Baseline snapshotted before review (branch review-gl.filter.pa, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.

test_that("SNP filter keeps exactly the private/fixed union", {
  pn <- popNames(testset.gl)
  f <- gl.filter.pa(testset.gl, pop1 = pn[1], pop2 = pn[2], verbose = 0)
  p1 <- as.matrix(testset.gl[pop(testset.gl) == pn[1], ])
  p2 <- as.matrix(testset.gl[pop(testset.gl) == pn[2], ])
  a1 <- colMeans(p1, na.rm = TRUE) / 2
  a2 <- colMeans(p2, na.rm = TRUE) / 2
  manual <- sum((a2 == 0 & a1 != 0) | (a2 == 1 & a1 != 1) |
                (a1 == 0 & a2 != 0) | (a1 == 1 & a2 != 1), na.rm = TRUE)
  expect_equal(nLoc(f), manual)
  expect_equal(nrow(f@other$loc.metrics), nLoc(f))
  # invers returns the complement
  fi <- gl.filter.pa(testset.gl, pop1 = pn[1], pop2 = pn[2],
                     invers = TRUE, verbose = 0)
  expect_equal(nLoc(fi), nLoc(testset.gl) - manual)
})

test_that("SilicoDArT arithmetic finds all private loci", {
  # [approved diff G1] baseline: allele frequencies were halved for
  # SilicoDArT too (max 0.5), so presence-fixed private alleles were
  # missed: 45 loci kept where the true union is 66 on testset.gs.
  pn <- popNames(testset.gs)
  f <- gl.filter.pa(testset.gs, pop1 = pn[1], pop2 = pn[2], verbose = 0)
  s1 <- as.matrix(testset.gs[pop(testset.gs) == pn[1], ])
  s2 <- as.matrix(testset.gs[pop(testset.gs) == pn[2], ])
  a1 <- colMeans(s1, na.rm = TRUE)
  a2 <- colMeans(s2, na.rm = TRUE)
  manual <- sum((a2 == 0 & a1 != 0) | (a2 == 1 & a1 != 1) |
                (a1 == 0 & a2 != 0) | (a1 == 1 & a2 != 1), na.rm = TRUE)
  expect_equal(nLoc(f), manual)  # [approved diff G1]
  expect_equal(manual, 66)
})

test_that("a bogus population name fails informatively", {
  # [approved diff G2] baseline: "'data' must be of a vector type,
  # was 'NULL'" from as.matrix(NULL).
  expect_error(
    gl.filter.pa(testset.gl, pop1 = "nonesuch",
                 pop2 = popNames(testset.gl)[2], verbose = 0),
    "nonesuch")  # [approved diff G2]
})

test_that("the filtered object returns invisibly with history appended", {
  pn <- popNames(testset.gl)
  # [approved diff G3] baseline: the object returned visibly.
  vf <- withVisible(gl.filter.pa(testset.gl, pop1 = pn[1], pop2 = pn[2],
                                 verbose = 0))
  expect_false(vf$visible)  # [approved diff G3]
  h <- vf$value@other$history
  expect_true(is.call(h[[length(h)]]))
})
