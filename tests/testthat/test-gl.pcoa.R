# Characterization tests for gl.pcoa
# Baseline snapshotted before review (review-gl.pcoa), upstream/dev ddaed27.
# Assertions marked [approved F<n>] were flipped in Phase C per the
# custodian's approval of findings 1, 3, 4, 5, 6 (2026-09-06).

test_that("glPca path: nfactors respected, eig full length, values stable", {
  pdf(NULL); on.exit(dev.off())
  for (nf in c(2, 5, 10)) {
    p <- gl.pcoa(testset2.gl, nfactors = nf, plot.out = FALSE, verbose = 0)
    expect_equal(dim(p$scores), c(274, nf))
    expect_equal(dim(p$loadings), c(611, nf))
    expect_length(p$eig, 273)
  }
  p5 <- gl.pcoa(testset2.gl, nfactors = 5, plot.out = FALSE, verbose = 0)
  expect_equal(round(p5$eig[1:3], 6), c(1.382214, 1.344892, 1.114741))
  expect_equal(round(sum(abs(p5$scores)), 3), 1186.296)
  expect_s3_class(p5, "glPca")
})

test_that("glPca path: verbose = 3 summary on testset2.gl has no NA lines", {
  pdf(NULL); on.exit(dev.off())
  out <- capture.output(
    invisible(gl.pcoa(testset2.gl, plot.out = FALSE, verbose = 3)))
  expect_true(any(grepl("27 informative dimensions", out)))
  expect_false(any(grepl("NA % of the total variance", out)))
})

test_that("low-rank object (<3 informative axes): verbose = 3 output", {
  pdf(NULL); on.exit(dev.off())
  lr <- testset2.gl[1:4, ]
  out <- capture.output(plr <- gl.pcoa(lr, plot.out = FALSE, verbose = 3))
  expect_length(plr$eig, 3)
  # [approved F5] axis-combination lines are now printed only when enough
  # informative axes exist; the NA line is gone and Axis 1-3 is suppressed.
  expect_false(any(grepl("NA % of the total variance", out)))
  expect_true(any(grepl("Axis 1 and 2 combined explain 61.7", out)))
  expect_false(any(grepl("Axis 1-3 combined", out)))
})

test_that("dist path: enough entities, nfactors honoured, values stable", {
  pdf(NULL); on.exit(dev.off())
  set.seed(42)
  m10 <- matrix(sample(0:2, 10 * 40, replace = TRUE), nrow = 10,
                dimnames = list(paste0("I", 1:10), NULL))
  D10 <- dist(m10)
  pd <- gl.pcoa(D10, nfactors = 5, plot.out = FALSE, verbose = 0)
  expect_equal(dim(pd$scores), c(10, 5))
  expect_length(pd$eig, 9)
  expect_equal(round(pd$eig[1:3], 4), c(44.7602, 38.2836, 33.5782))
  # correction = "none": ape::pcoa returns no vectors.cor, loadings are NULL
  expect_null(pd$loadings)
})

test_that("dist path: fewer entities than nfactors + 1", {
  pdf(NULL); on.exit(dev.off())
  set.seed(42)
  m4 <- matrix(sample(0:2, 4 * 40, replace = TRUE), nrow = 4,
               dimnames = list(paste0("I", 1:4), NULL))
  D4 <- dist(m4)
  # [approved F6] previously died with "subscript out of bounds"; now clamps
  # nfactors to the available axes and warns.
  out <- capture.output(pd4 <- gl.pcoa(D4, nfactors = 5, plot.out = FALSE,
                                       verbose = 2))
  expect_true(any(grepl("Warning.*nfactors", out)))
  expect_equal(dim(pd4$scores), c(4, 3))
  expect_length(pd4$eig, 3)
})

test_that("dist path: verbose = 3 with <3 positive eigenvalues", {
  pdf(NULL); on.exit(dev.off())
  set.seed(7)
  m3 <- matrix(sample(0:2, 3 * 40, replace = TRUE), nrow = 3,
               dimnames = list(paste0("I", 1:3), NULL))
  D3 <- dist(m3)
  out <- capture.output(p3 <- gl.pcoa(D3, nfactors = 2, plot.out = FALSE,
                                      verbose = 3))
  expect_equal(dim(p3$scores), c(3, 2))
  # [approved F5] the Axis 1-3 line is suppressed when only two positive
  # eigenvalues exist (guard was length(eig.top >= 3), a precedence bug).
  expect_false(any(grepl("NA % of the total variance", out)))
  expect_true(any(grepl("Axis 1 and 2 combined explain 100", out)))
  expect_false(any(grepl("Axis 1-3 combined", out)))
})

test_that("FBM path: nfactors and eig length", {
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off())
  fb <- gl.gen2fbm(testset2.gl)
  pf <- gl.pcoa(fb, nfactors = 5, plot.out = FALSE, verbose = 0)
  # [approved F1] nfactors was silently ignored on the big_SVD path (scores
  # came back 274 x 273); scores/loadings are now truncated to nfactors,
  # $eig stays full length as on the glPca path.
  expect_equal(dim(pf$scores), c(274, 5))
  expect_equal(dim(pf$loadings), c(611, 5))
  expect_length(pf$eig, 273)
  expect_equal(round(pf$eig[1:3], 6), c(2.259724, 2.184270, 1.652706))
})
