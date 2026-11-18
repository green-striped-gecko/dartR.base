# Characterization tests for utils.dist.ind.snp
# Baseline snapshotted before review (distance-wave review, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.

make_g <- function() {
  set.seed(21)
  m <- matrix(sample(0:2, 10 * 30, replace = TRUE), nrow = 10)
  m[sample(length(m), 20)] <- NA
  g <- new("genlight", gen = m, ploidy = 2)
  indNames(g) <- paste0("i", 1:10); locNames(g) <- paste0("L", 1:30)
  pop(g) <- factor(rep(c("A", "B"), each = 5))
  suppressWarnings(gl.compliance.check(g, verbose = 0))
}
relabel <- function(g) {
  m <- 2 - as.matrix(g)
  g2 <- new("genlight", gen = m, ploidy = 2)
  indNames(g2) <- indNames(g); locNames(g2) <- locNames(g); pop(g2) <- pop(g)
  suppressWarnings(gl.compliance.check(g2, verbose = 0))
}

test_that("euclidean and manhattan match manual recomputations", {
  g <- make_g(); mm <- as.matrix(g)
  me <- matrix(0, 10, 10); mman <- matrix(0, 10, 10)
  for (i in 1:9) for (j in (i + 1):10) {
    d <- mm[i, ] - mm[j, ]; ok <- !is.na(d)
    me[i, j] <- me[j, i] <- sqrt(sum(d[ok]^2))
    mman[i, j] <- mman[j, i] <- sum(abs(d[ok])) / (2 * sum(ok))
  }
  o <- capture.output({
    eu <- as.matrix(utils.dist.ind.snp(g, method = "Euclidean", verbose = 0))
    ma <- as.matrix(utils.dist.ind.snp(g, method = "Manhattan", verbose = 0))
  })
  expect_equal(unname(eu), unname(me))
  expect_equal(unname(ma), unname(mman))
})

test_that("simple and absolute are invariant under allele relabelling", {
  # [approved diff D1] baseline: both methods scored both-hom-REF pairs as
  # sharing no alleles while both-hom-ALT shared both, so recoding which
  # allele is reference changed the distances (max discrepancy 0.167 on
  # this fixture).
  g <- make_g(); g2 <- relabel(g)
  o <- capture.output({
    s1 <- as.matrix(utils.dist.ind.snp(g, method = "Simple", verbose = 0))
    s2 <- as.matrix(utils.dist.ind.snp(g2, method = "Simple", verbose = 0))
    a1 <- as.matrix(utils.dist.ind.snp(g, method = "Absolute", verbose = 0))
    a2 <- as.matrix(utils.dist.ind.snp(g2, method = "Absolute", verbose = 0))
  })
  expect_equal(s1, s2)  # [approved diff D1]
  expect_equal(a1, a2)  # [approved diff D1]
})

test_that("verbose 0 is silent and messages print once", {
  g <- make_g()
  # [approved diff D2] baseline: an ungated progress cat leaked 1 line at
  # verbose 0, and the duplicated tail printed Completed/Returning twice.
  o0 <- capture.output(d <- utils.dist.ind.snp(g, verbose = 0))
  expect_equal(length(o0), 0)  # [approved diff D2]
  o2 <- capture.output(d <- utils.dist.ind.snp(g, verbose = 2))
  expect_equal(sum(grepl("Completed", o2)), 1)  # [approved diff D2]
})
