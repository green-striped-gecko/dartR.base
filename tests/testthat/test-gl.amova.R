# Characterization tests for gl.amova — Phase A baseline (pre-change at
# ddaed27), updated in Phase C: assertions tagged [approved Fn] flipped
# with the matching approved finding. See
# function-review/reports/dartR.base/gl.amova.md for findings.

test_that("gl.amova pins components, phi and p on testset.gl (seeded)", {
  set.seed(42)
  res <- gl.amova(testset.gl, permutations = 99, verbose = 0)

  # structure
  expect_s3_class(res, "amova")
  expect_named(res, c("tab", "varcoef", "varcomp", "call"))
  expect_equal(rownames(res$tab), c("pop.names", "Error", "Total"))

  # pinned table and components
  expect_equal(res$tab$df, c(29, 220, 249))
  expect_equal(res$tab$SSD,
               c(0.111257489413457, 0.00505738501540252, 0.11631487442886),
               tolerance = 1e-10)
  expect_equal(res$varcomp$sigma2,
               c(0.000459523792989282, 2.29881137063751e-05),
               tolerance = 1e-10)
  # pegas p = sum(rand >= obs)/(nperm + 1): zero attainable
  expect_equal(res$varcomp$P.value, c(0, NA))

  # derived phi_ST
  phi <- res$varcomp$sigma2[1] / sum(res$varcomp$sigma2)
  expect_equal(phi, 0.952357, tolerance = 1e-6)
})

test_that("gl.amova is deterministic under set.seed and honours permutations", {
  set.seed(42)
  r1 <- gl.amova(testset.gl, permutations = 99, verbose = 0)
  set.seed(42)
  r2 <- gl.amova(testset.gl, permutations = 99, verbose = 0)
  expect_identical(r1$varcomp, r2$varcomp)

  # sigma2 does not depend on the permutation count
  set.seed(1)
  r3 <- gl.amova(testset.gl, permutations = 9, verbose = 0)
  expect_equal(r3$varcomp$sigma2, r1$varcomp$sigma2, tolerance = 1e-12)
})

test_that("gl.amova equals direct pegas::amova on stamppNeisD (delegation)", {
  x <- testset.gl
  class(x) <- "genlight"
  dd <- StAMPP::stamppNeisD(x, FALSE)
  expect_identical(rownames(dd), indNames(testset.gl))
  pops <- factor(as.character(pop(testset.gl)))
  ddist <- stats::as.dist(dd)
  set.seed(42)
  ind <- pegas::amova(ddist ~ pops, nperm = 99)
  set.seed(42)
  res <- gl.amova(testset.gl, permutations = 99, verbose = 0)
  expect_equal(res$tab$SSD, ind$tab$SSD, tolerance = 1e-12)
  expect_equal(res$varcomp$sigma2, ind$varcomp$sigma2, tolerance = 1e-12)
})

test_that("gl.amova matches hand-computed one-level AMOVA (Excoffier 1992)", {
  x <- testset.gl
  class(x) <- "genlight"
  dd <- StAMPP::stamppNeisD(x, FALSE)
  pops <- factor(as.character(pop(testset.gl)))
  D2 <- as.matrix(stats::as.dist(dd))^2
  N <- nrow(D2)
  g <- nlevels(pops)
  ssd_tot <- sum(D2[upper.tri(D2)]) / N
  ssd_wit <- 0
  for (p in levels(pops)) {
    idx <- which(pops == p)
    sub <- D2[idx, idx, drop = FALSE]
    ssd_wit <- ssd_wit + sum(sub[upper.tri(sub)]) / length(idx)
  }
  ssd_amg <- ssd_tot - ssd_wit
  msd_wit <- ssd_wit / (N - g)
  msd_amg <- ssd_amg / (g - 1)
  n0 <- (N - sum(tapply(seq_along(pops), pops, length)^2) / N) / (g - 1)
  sig2_w <- msd_wit
  sig2_a <- (msd_amg - msd_wit) / n0

  set.seed(42)
  res <- gl.amova(testset.gl, permutations = 9, verbose = 0)
  expect_equal(res$tab$SSD, c(ssd_amg, ssd_wit, ssd_tot), tolerance = 1e-9)
  expect_equal(res$varcomp$sigma2, c(sig2_a, sig2_w), tolerance = 1e-9)
})

test_that("gl.amova accepts a supplied distance (matrix and dist) identically", {
  x <- testset.gl
  class(x) <- "genlight"
  dd <- StAMPP::stamppNeisD(x, FALSE)
  set.seed(42)
  internal <- gl.amova(testset.gl, permutations = 9, verbose = 0)
  set.seed(42)
  supplied_m <- gl.amova(testset.gl, distance = dd, permutations = 9,
                         verbose = 0)
  set.seed(42)
  supplied_d <- gl.amova(testset.gl, distance = stats::as.dist(dd),
                         permutations = 9, verbose = 0)
  expect_identical(internal$varcomp$sigma2, supplied_m$varcomp$sigma2)
  expect_equal(supplied_d$varcomp$sigma2, internal$varcomp$sigma2,
               tolerance = 1e-12)
})

test_that("supplied distance is aligned to indNames(x) by its labels", {
  # [approved F1] flipped: a row-permuted matrix with correct labels was
  # previously used positionally (silently wrong components); it is now
  # aligned by name and gives the same results as the ordered matrix
  x <- testset.gl
  class(x) <- "genlight"
  dd <- StAMPP::stamppNeisD(x, FALSE)
  set.seed(7)
  perm <- sample(nrow(dd))
  ddp <- dd[perm, perm]
  rownames(ddp) <- rownames(dd)[perm]  # labels are correct, order is not
  colnames(ddp) <- rownames(dd)[perm]
  set.seed(42)
  good <- gl.amova(testset.gl, permutations = 9, verbose = 0)
  set.seed(42)
  aligned <- gl.amova(testset.gl, distance = ddp, permutations = 9,
                      verbose = 0)
  expect_equal(aligned$varcomp$sigma2, good$varcomp$sigma2,
               tolerance = 1e-12)
})

test_that("unlabelled and mislabelled distances are fatal", {
  # [approved F1] a distance without individual names cannot be matched
  # to the genlight object and now errors informatively
  x <- testset.gl
  class(x) <- "genlight"
  dd <- StAMPP::stamppNeisD(x, FALSE)
  ddu <- unname(dd)
  expect_error(
    gl.amova(testset.gl, distance = ddu, permutations = 9, verbose = 0),
    "no individual names"
  )
  ddw <- dd
  rownames(ddw)[1] <- "not_an_individual"
  colnames(ddw) <- rownames(ddw)
  expect_error(
    gl.amova(testset.gl, distance = ddw, permutations = 9, verbose = 0),
    "do not match indNames"
  )
})

test_that("wrong-dimension distance fails informatively", {
  # [approved F1] flipped: previously died inside pegas with bare
  # "subscript out of bounds"
  x <- testset.gl
  class(x) <- "genlight"
  dd <- StAMPP::stamppNeisD(x, FALSE)
  expect_error(
    gl.amova(testset.gl, distance = dd[1:10, 1:10], permutations = 9,
             verbose = 0),
    "covers 10 individuals"
  )
})

test_that("SilicoDArT input is rejected (SNP-only contract)", {
  # [approved F4] flipped: testset.gs previously ran and returned
  # components on an unvalidated presence-absence distance basis
  expect_error(
    gl.amova(testset.gs, permutations = 9, verbose = 0),
    "works only with SNP data"
  )
})

test_that("gl.amova two-population subset pins", {
  x2 <- gl.keep.pop(testset.gl,
                    pop.list = c("EmmacBurnBara", "EmmacMaclGeor"),
                    verbose = 0)
  set.seed(42)
  r2 <- gl.amova(x2, permutations = 99, verbose = 0)
  expect_equal(r2$tab$df[1], 1)
  expect_equal(r2$varcomp$sigma2,
               c(6.28029064921488e-05, 1.926336395e-05),
               tolerance = 1e-10)
  expect_equal(r2$varcomp$P.value, c(0, NA))
})

test_that("single-population input is fatal", {
  # [approved F3] flipped: previously ran to completion and returned an
  # all-NaN table silently
  x1 <- gl.keep.pop(testset.gl, pop.list = "EmmacBurnBara", verbose = 0)
  expect_error(
    gl.amova(x1, permutations = 9, verbose = 0),
    "at least two populations"
  )
})

test_that("non-finite distances are fatal and name the individuals", {
  # [approved F6] flipped: an all-NA individual previously propagated to
  # an all-NaN table with only pegas's generic warning
  x <- testset.gl
  class(x) <- "genlight"
  dd <- StAMPP::stamppNeisD(x, FALSE)
  colnames(dd) <- rownames(dd)
  dd[2, 1] <- dd[1, 2] <- NaN
  expect_error(
    gl.amova(testset.gl, distance = dd, permutations = 9, verbose = 0),
    rownames(dd)[1]
  )
})

test_that("gl.amova components are invariant to individual order and
          pop label alphabet", {
  x4 <- testset.gl
  pop(x4) <- factor(paste0(
    ifelse(as.integer(pop(x4)) %% 2 == 0, "zzz_", "aaa_"),
    as.character(pop(x4))))
  set.seed(42)
  rL <- gl.amova(x4, permutations = 9, verbose = 0)
  set.seed(11)
  ord <- sample(nInd(x4))
  x5 <- x4[ord, ]
  set.seed(42)
  rS <- gl.amova(x5, permutations = 9, verbose = 0)
  expect_equal(rL$varcomp$sigma2, rS$varcomp$sigma2, tolerance = 1e-9)
})

test_that("gl.amova is silent at verbose = 0 and summarises at verbose = 3", {
  set.seed(1)
  out <- capture.output(res <- gl.amova(testset.gl, permutations = 9,
                                        verbose = 0))
  expect_length(out, 0)
  expect_s3_class(res, "amova")

  # [approved F7] a results summary now prints at verbose >= 3
  set.seed(1)
  out3 <- capture.output(res3 <- gl.amova(testset.gl, permutations = 9,
                                          verbose = 3))
  expect_true(any(grepl("Phi_ST", out3)))
  expect_true(any(grepl("Variance components", out3)))
})
