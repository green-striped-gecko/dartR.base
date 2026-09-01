# Characterization tests for utils.basic.stats
# Baseline snapshotted before review (kernel-wave review, dev at ddaed27).
# Assertions marked [approved diff] were flipped in Phase C.

make_g <- function(perpop_allna = FALSE) {
  set.seed(11)
  m <- matrix(sample(0:2, 30 * 40, replace = TRUE), nrow = 30)
  if (perpop_allna) {
    m[1:10, 3] <- NA
    m[11:20, c(7, 8)] <- NA
  }
  g <- new("genlight", gen = m, ploidy = 2)
  indNames(g) <- paste0("i", 1:30); locNames(g) <- paste0("L", 1:40)
  pop(g) <- factor(rep(c("A", "B", "C"), each = 10))
  suppressWarnings(gl.compliance.check(g, verbose = 0))
}

test_that("matches hierfstat::basic.stats exactly on complete data", {
  skip_if_not_installed("hierfstat")
  g <- make_g()
  ub <- utils.basic.stats(g)
  hb <- hierfstat::basic.stats(hierfstat::genind2hierfstat(gl2gi(g, verbose = 0)))
  expect_equal(unname(ub$overall[1:9]), unname(round(hb$overall[1:9], 4)),
               tolerance = 1e-3)
  expect_equal(unname(ub$perloc$Hs), unname(round(hb$perloc$Hs, 4)),
               tolerance = 1e-3)
})

test_that("loci absent from one population match hierfstat", {
  skip_if_not_installed("hierfstat")
  # [approved diff K1] baseline: a locus with zero data in any population
  # made the harmonic mean sample size 0, so Hs/Ht/Fst went NaN for that
  # locus and it silently dropped from the overall averages (hierfstat
  # computes it from the populations that have data). The doc claims
  # exact equivalence.
  g <- make_g(perpop_allna = TRUE)
  ub <- utils.basic.stats(g)
  hb <- hierfstat::basic.stats(hierfstat::genind2hierfstat(gl2gi(g, verbose = 0)))
  expect_equal(ub$perloc$Hs[3], round(hb$perloc$Hs[3], 4), tolerance = 1e-3)  # [approved diff K1]
  expect_equal(unname(ub$overall[1:9]),
    unname(round(hb$overall[1:9], 4)), tolerance = 1e-3)  # [approved diff K1]
})

test_that("a single-individual population fails with a message", {
  g <- make_g()
  gg <- g[1:11, ]
  pop(gg) <- factor(c(rep("A", 10), "B"))
  e <- tryCatch(utils.basic.stats(gg), error = function(e) conditionMessage(e))
  # [approved diff K2] baseline: cat(error()) + bare stop() left the error
  # message empty.
  expect_true(grepl("individual", e))  # [approved diff K2]
})
