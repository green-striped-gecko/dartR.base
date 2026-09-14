# Tests for gl.sim.crosses (dartR.base)
#
# Captured during the function-review campaign at upstream/dev = ddaed27 as a
# characterization baseline, then updated in the same change that applied the
# approved findings. Every expectation that changed carries an "[approved Fn]"
# marker naming the finding that changed it; everything else is unchanged from
# the baseline and must keep passing.
#
# Report: function-review/reports/dartR.base/gl.sim.crosses.md

# ---- fixtures ---------------------------------------------------------------

# Minimal genlight built straight from a dosage matrix. Deliberately not passed
# through gl.compliance.check: the engineered fixtures below are one- and
# two-locus objects that compliance checking rejects, and the function reads
# only as.matrix(), locNames(), indNames() and nInd() from its inputs.
sc_fixture <- function(mat, prefix = "P") {
  rownames(mat) <- paste0(prefix, seq_len(nrow(mat)))
  colnames(mat) <- paste0("L", seq_len(ncol(mat)))
  new("genlight",
      gen = mat,
      ind.names = rownames(mat),
      loc.names = colnames(mat),
      ploidy = rep(2, nrow(mat)))
}

sc_run <- function(...) {
  invisible(capture.output(res <- suppressWarnings(gl.sim.crosses(...))))
  res
}

# ---- documented parameters --------------------------------------------------

test_that("the documented `n` parameter retains n offspring", {
  # [approved F2] was a DEFECT PIN: n was read only by a warning test and the
  # full brood was always returned (15 here).
  mums <- sc_fixture(matrix(c(1, 1, 0, 2,
                              1, 0, 1, 2,
                              0, 1, 1, 2), nrow = 3, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(1, 1, 2, 0,
                              0, 1, 1, 0,
                              1, 1, 0, 0), nrow = 3, byrow = TRUE), "F")
  set.seed(1)
  out <- sc_run(dads, mums, broodsize = 5, n = 2, verbose = 0,
                compliance.check = FALSE)
  expect_equal(nInd(out), 2)
  expect_equal(indNames(out), c("Po_1", "Po_2"))
  set.seed(1)
  out6 <- sc_run(dads, mums, broodsize = 5, n = 6, verbose = 0,
                 compliance.check = FALSE)
  expect_equal(nInd(out6), 6)
  # genotypes, names and metadata are subset together
  expect_equal(nrow(out6@other$ind.metrics), 6)
  expect_equal(nrow(as.matrix(out6)), 6)
})

test_that("[approved F2] the default n resolves to the lesser of 1000 and the brood total", {
  mums <- sc_fixture(matrix(c(1, 1, 0, 2,
                              1, 0, 1, 2,
                              0, 1, 1, 2), nrow = 3, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(1, 1, 2, 0,
                              0, 1, 1, 0,
                              1, 1, 0, 0), nrow = 3, byrow = TRUE), "F")
  out <- sc_run(dads, mums, broodsize = 5, verbose = 0,
                compliance.check = FALSE)
  expect_equal(nInd(out), 15)                 # brood total is the lesser
  big <- sc_run(dads, mums, broodsize = 400, verbose = 0,
                compliance.check = FALSE)
  expect_equal(nInd(big), 1000)               # 1000 is the lesser
})

test_that("[approved F2] a non-positive or fractional n is fatal", {
  mums <- sc_fixture(matrix(c(2, 1, 2, 1), nrow = 2, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(2, 1, 0, 1), nrow = 2, byrow = TRUE), "F")
  expect_error(gl.sim.crosses(dads, mums, broodsize = 2, n = 0, verbose = 0),
               "positive integer")
  expect_error(gl.sim.crosses(dads, mums, broodsize = 2, n = 2.5, verbose = 0),
               "positive integer")
})

test_that("error.check = FALSE completes", {
  # [approved F3] was a DEFECT PIN: `noff` was computed inside the error-check
  # block and used by the constructor, so the setting @details recommends for
  # simulations failed with "object 'noff' not found".
  mums <- sc_fixture(matrix(c(2, 1, 2, 1), nrow = 2, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(2, 1, 0, 1), nrow = 2, byrow = TRUE), "F")
  out <- sc_run(dads, mums, broodsize = 5, error.check = FALSE, verbose = 0)
  expect_equal(nInd(out), 10)
  expect_equal(indNames(out), paste0("Po_", 1:10))
  expect_true(all(nzchar(indNames(out))))
})

test_that("brood size determines the number of offspring returned", {
  mums <- sc_fixture(matrix(c(2, 1, 2, 1), nrow = 2, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(2, 1, 0, 1), nrow = 2, byrow = TRUE), "F")
  for (bs in c(1, 2, 7)) {
    out <- sc_run(dads, mums, broodsize = bs, verbose = 0,
                  compliance.check = FALSE)
    expect_equal(nInd(out), 2 * bs)
    expect_equal(indNames(out), paste0("Po_", seq_len(2 * bs)))
  }
})

test_that("non-positive or fractional broodsize falls back to 10", {
  # [approved F6] was a DEFECT PIN: `for (i in 1:broodsize)` counted DOWN when
  # broodsize <= 0, `noff` then disagreed with nrow(offmat), and every
  # individual name came back an empty string. The documented "Set to 10" is
  # now assigned rather than merely announced.
  mums <- sc_fixture(matrix(c(2, 1, 2, 1), nrow = 2, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(2, 1, 0, 1), nrow = 2, byrow = TRUE), "F")
  for (bs in c(0, -1, 2.7)) {
    out <- sc_run(dads, mums, broodsize = bs, verbose = 0,
                  compliance.check = FALSE)
    expect_equal(nInd(out), 20)                    # 2 mothers x 10
    expect_true(all(nzchar(indNames(out))))
    expect_equal(indNames(out), paste0("Po_", 1:20))
  }
  txt <- capture.output(suppressWarnings(
    res <- gl.sim.crosses(dads, mums, broodsize = 0, verbose = 1,
                   compliance.check = FALSE)))
  expect_true(any(grepl("Brood size must be a positive integer", txt)))
})

test_that("out-of-range sexratio falls back to 0.5", {
  # [approved F12] was a DEFECT PIN: the guard printed "Set to 0.5" but never
  # assigned it, so sexratio = 1.7 returned an all-female, single-level cohort.
  mums <- sc_fixture(matrix(c(2, 1, 2, 1), nrow = 2, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(2, 1, 0, 1), nrow = 2, byrow = TRUE), "F")
  set.seed(4)
  out <- sc_run(dads, mums, broodsize = 50, sexratio = 1.7, verbose = 0,
                compliance.check = FALSE)
  expect_equal(levels(out@other$ind.metrics$sex), c("female", "male"))
  expect_equal(mean(out@other$ind.metrics$sex == "female"), 0.5,
               tolerance = 0.15)
  txt <- capture.output(suppressWarnings(
    res <- gl.sim.crosses(dads, mums, broodsize = 2, sexratio = 1.7, verbose = 1,
                   compliance.check = FALSE)))
  expect_true(any(grepl("Sex ratio must be in the range 0 to 1", txt)))
})

test_that("sexratio is the expected proportion of females", {
  mums <- sc_fixture(matrix(c(2, 1, 2, 1), nrow = 2, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(2, 1, 0, 1), nrow = 2, byrow = TRUE), "F")
  set.seed(4)
  out <- sc_run(dads, mums, broodsize = 500, sexratio = 0.9, verbose = 0,
                compliance.check = FALSE)
  sx <- table(out@other$ind.metrics$sex)
  expect_equal(as.numeric(sx["female"] / sum(sx)), 0.9, tolerance = 0.05)
  # [approved F12] both levels are declared even for a single-sexed cohort
  out1 <- sc_run(dads, mums, broodsize = 20, sexratio = 1, verbose = 0,
                 compliance.check = FALSE)
  expect_equal(levels(out1@other$ind.metrics$sex), c("female", "male"))
})

# ---- Mendelian correctness --------------------------------------------------

test_that("offspring dosages follow Mendelian expectations at the margin", {
  # One pair, six loci covering every cross type.
  #  L1 0x0 -> all 0     L2 2x2 -> all 2     L3 0x2 -> all 1
  #  L4 1x1 -> 1:2:1     L5 0x1 -> 1:1 (0/1) L6 1x2 -> 1:1 (1/2)
  # [approved F2] n = 4000 is now passed explicitly; the default would cap the
  # brood at 1000. The marginal ratios themselves are unchanged.
  mums <- sc_fixture(matrix(c(0, 2, 0, 1, 0, 1), nrow = 1), "M")
  dads <- sc_fixture(matrix(c(0, 2, 2, 1, 1, 2), nrow = 1), "F")
  set.seed(101)
  out <- sc_run(dads, mums, broodsize = 4000, n = 4000, verbose = 0,
                compliance.check = FALSE)
  m <- as.matrix(out)
  expect_equal(nrow(m), 4000)
  expect_true(all(m[, 1] == 0))
  expect_true(all(m[, 2] == 2))
  expect_true(all(m[, 3] == 1))
  p4 <- prop.table(table(factor(m[, 4], levels = 0:2)))
  expect_equal(as.numeric(p4), c(0.25, 0.5, 0.25), tolerance = 0.03)
  p5 <- prop.table(table(factor(m[, 5], levels = 0:2)))
  expect_equal(as.numeric(p5), c(0.5, 0.5, 0), tolerance = 0.03)
  p6 <- prop.table(table(factor(m[, 6], levels = 0:2)))
  expect_equal(as.numeric(p6), c(0, 0.5, 0.5), tolerance = 0.03)
  # only 0/1/2 are produced
  expect_true(all(m %in% c(0, 1, 2)))
})

test_that("offspring of monomorphic parents are fixed, and single-locus objects work", {
  mums <- sc_fixture(matrix(c(0, 2, 0, 2), nrow = 2, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(0, 2, 0, 2), nrow = 2, byrow = TRUE), "F")
  out <- sc_run(dads, mums, broodsize = 3, verbose = 0, compliance.check = FALSE)
  m <- as.matrix(out)
  expect_true(all(m[, 1] == 0))
  expect_true(all(m[, 2] == 2))
  m1 <- sc_fixture(matrix(c(1, 1), nrow = 2), "M")
  d1 <- sc_fixture(matrix(c(1, 1), nrow = 2), "F")
  set.seed(5)
  o1 <- sc_run(d1, m1, broodsize = 3, verbose = 0, compliance.check = FALSE)
  expect_equal(nLoc(o1), 1)
  expect_equal(nInd(o1), 6)
  expect_true(all(as.matrix(o1) %in% c(0, 1, 2)))
})

# ---- gamete independence (the sampling defect) ------------------------------

test_that("het loci within one parent segregate independently", {
  # [approved F1] was a DEFECT PIN. The old idiom
  # `ifelse(mmat == 1, sample(c(0, 2), mhet, replace = TRUE), mmat)` drew only
  # `mhet` values and let ifelse() recycle them over the whole matrix, so het
  # calls whose column-major positions were congruent mod mhet always
  # transmitted the same allele. Here the mother is het at L1 and L3 with
  # mhet = 2, so both took draw 1 and L1 == L3 in 500/500 offspring.
  mums <- sc_fixture(matrix(c(1, 0, 1, 2), nrow = 1), "M")
  dads <- sc_fixture(matrix(c(0, 0, 0, 2), nrow = 1), "F")
  set.seed(7)
  out <- sc_run(dads, mums, broodsize = 500, n = 500, verbose = 0,
                compliance.check = FALSE)
  m <- as.matrix(out)
  expect_false(all(m[, 1] == m[, 3]))
  expect_equal(mean(m[, 1] == m[, 3]), 0.5, tolerance = 0.1)
  # the marginal ratio at each locus is still 1:1
  expect_equal(mean(m[, 1] == 1), 0.5, tolerance = 0.1)
  expect_equal(mean(m[, 3] == 1), 0.5, tolerance = 0.1)
})

test_that("unrelated parents draw gametes independently", {
  # [approved F1] was a DEFECT PIN. Three mothers, two loci. M1 and M3 are het
  # at L1 and nothing else is het, so mhet = 2 and their two het cells
  # (column-major positions 1 and 3) both took draw 1: agreement was 300/300.
  mums <- sc_fixture(matrix(c(1, 2,
                              0, 2,
                              1, 2), nrow = 3, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(0, 2,
                              0, 2,
                              0, 2), nrow = 3, byrow = TRUE), "F")
  set.seed(9)
  out <- sc_run(dads, mums, broodsize = 300, n = 900, verbose = 0,
                compliance.check = FALSE)
  m <- as.matrix(out)
  r1 <- seq(1, nrow(m), by = 3)   # offspring of mother 1
  r3 <- seq(3, nrow(m), by = 3)   # offspring of mother 3
  expect_false(all(m[r1, 1] == m[r3, 1]))
  expect_equal(mean(m[r1, 1] == m[r3, 1]), 0.5, tolerance = 0.1)
})

test_that("every heterozygous call consumes its own random draw", {
  # [approved F1] Structural pin on testset.gl: 10 mothers x 255 loci carry 29
  # het calls. Under the old ifelse() recycling those 29 cells addressed only
  # 14 distinct draws; one draw per het cell gives 29.
  out1 <- capture.output(fems <- gl.keep.pop(testset.gl, pop.list = "Female",
                                            as.pop = "sex", verbose = 0))
  set.seed(32)
  mo <- suppressWarnings({
    o <- capture.output(v <- gl.keep.ind(fems, sample(indNames(fems), 10),
                                        verbose = 0)); v })
  mmat <- as.matrix(mo)
  mhet <- sum(mmat == 1, na.rm = TRUE)
  k <- which(mmat == 1)
  old.index <- ((k - 1) %% mhet) + 1     # the draw each het cell used to get
  new.index <- seq_along(k)              # the draw each het cell now gets
  expect_equal(mhet, 29)
  expect_equal(length(unique(old.index)), 14)
  expect_equal(length(unique(new.index)), mhet)

  # and the two cells the old code forced to share a draw now disagree about
  # half the time. Fathers are homozygous reference, so the offspring dosage
  # reads the maternal allele directly.
  dup <- old.index[duplicated(old.index)][1]
  pair <- k[old.index == dup][1:2]
  rc <- arrayInd(pair, dim(mmat))
  fa0 <- sc_fixture(matrix(0, nrow = nInd(mo), ncol = nLoc(mo)), "Fa")
  locNames(fa0) <- locNames(mo)
  set.seed(77)
  off <- sc_run(fa0, mo, broodsize = 200, n = 200 * nInd(mo), verbose = 0,
                compliance.check = FALSE)
  om <- as.matrix(off)
  a <- om[seq(rc[1, 1], nrow(om), by = nInd(mo)), rc[1, 2]]
  b <- om[seq(rc[2, 1], nrow(om), by = nInd(mo)), rc[2, 2]]
  expect_false(all(a == b))
  expect_equal(mean(a == b), 0.5, tolerance = 0.15)
})

test_that("full sibs are not identical: brood replicates draw independently", {
  mums <- sc_fixture(matrix(c(1, 1, 1, 1), nrow = 1), "M")
  dads <- sc_fixture(matrix(c(1, 1, 1, 1), nrow = 1), "F")
  set.seed(17)
  out <- sc_run(dads, mums, broodsize = 200, verbose = 0,
                compliance.check = FALSE)
  m <- as.matrix(out)
  expect_gt(length(unique(apply(m, 1, paste, collapse = "-"))), 1)
})

# ---- missing data -----------------------------------------------------------

test_that("a missing call in either parent yields a missing call in the offspring", {
  mums <- sc_fixture(matrix(c(NA, 1, NA, 1, 2), nrow = 1), "M")
  dads <- sc_fixture(matrix(c(1, NA, NA, 1, 0), nrow = 1), "F")
  set.seed(11)
  out <- sc_run(dads, mums, broodsize = 50, verbose = 0,
                compliance.check = FALSE)
  m <- as.matrix(out)
  expect_true(all(is.na(m[, 1])))    # mother NA
  expect_true(all(is.na(m[, 2])))    # father NA
  expect_true(all(is.na(m[, 3])))    # both NA
  expect_false(any(is.na(m[, 4])))
  expect_true(all(m[, 5] == 1))
})

test_that("the offspring NA pattern reproduces the parental NA pattern exactly", {
  set.seed(3)
  g <- matrix(sample(c(0, 1, 2, NA), 20 * 30, replace = TRUE,
                     prob = c(0.3, 0.3, 0.3, 0.1)), nrow = 20)
  mums <- sc_fixture(g, "M")
  dads <- sc_fixture(g, "F")
  set.seed(13)
  out <- sc_run(dads, mums, broodsize = 3, verbose = 0,
                compliance.check = FALSE)
  m <- as.matrix(out)
  pidx <- ((seq_len(nrow(m)) - 1) %% 20) + 1
  pmat <- as.matrix(mums)
  expect_equal(unname(is.na(m)), unname(is.na(pmat[pidx, , drop = FALSE])))
  # homozygous parents can only produce that homozygote
  hom0 <- !is.na(pmat[pidx, ]) & pmat[pidx, ] == 0
  expect_true(all(m[hom0] == 0))
  hom2 <- !is.na(pmat[pidx, ]) & pmat[pidx, ] == 2
  expect_true(all(m[hom2] == 2))
})

test_that("repeated parents do not corrupt the NA pattern", {
  # SNPbin duplicate-index negative control: the function reads its inputs
  # through as.matrix() and never invokes the SNPbin `[` subsetter.
  set.seed(3)
  g <- matrix(sample(c(0, 1, 2, NA), 6 * 30, replace = TRUE,
                     prob = c(0.3, 0.3, 0.3, 0.1)), nrow = 6)
  mums <- sc_fixture(rbind(g[1, , drop = FALSE], g[1, , drop = FALSE],
                           g[2, , drop = FALSE]), "M")
  dads <- sc_fixture(rbind(g[3, , drop = FALSE], g[3, , drop = FALSE],
                           g[4, , drop = FALSE]), "F")
  set.seed(99)
  out <- sc_run(dads, mums, broodsize = 4, verbose = 0,
                compliance.check = FALSE)
  m <- as.matrix(out)
  pidx <- ((seq_len(nrow(m)) - 1) %% 3) + 1
  expected <- is.na(as.matrix(mums)[pidx, , drop = FALSE]) |
    is.na(as.matrix(dads)[pidx, , drop = FALSE])
  expect_equal(unname(is.na(m)), unname(expected))
})

test_that("an all-NA locus survives as an all-NA locus", {
  mums <- sc_fixture(matrix(c(NA, 1, NA, 1), nrow = 2, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(NA, 1, NA, 1), nrow = 2, byrow = TRUE), "F")
  set.seed(19)
  out <- sc_run(dads, mums, broodsize = 2, verbose = 0,
                compliance.check = FALSE)
  m <- as.matrix(out)
  expect_true(all(is.na(m[, 1])))
  expect_equal(nLoc(out), 2)
})

# ---- pairing and cohort validation -----------------------------------------

test_that("mothers are paired with fathers by position, identically in every brood", {
  mums <- sc_fixture(matrix(c(2, 1,
                              2, 1), nrow = 2, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(2, 1,
                              0, 1), nrow = 2, byrow = TRUE), "F")
  set.seed(2)
  out <- sc_run(dads, mums, broodsize = 3, verbose = 0,
                compliance.check = FALSE)
  m <- as.matrix(out)
  expect_true(all(m[c(1, 3, 5), 1] == 2))   # M1 x F1, both hom-alt
  expect_true(all(m[c(2, 4, 6), 1] == 1))   # M2 x F2, hom-alt x hom-ref
  # [approved F11] the pairing is now recorded
  expect_equal(out@other$ind.metrics$mother, rep(c("M1", "M2"), 3))
  expect_equal(out@other$ind.metrics$father, rep(c("F1", "F2"), 3))
})

test_that("unequal cohorts fail with an informative error", {
  # [approved F7] was a DEFECT PIN: the base R message "non-conformable arrays"
  # gave no indication of which argument was wrong.
  m5 <- sc_fixture(matrix(rep(c(0, 1, 2, 1, 0), 5), nrow = 5, byrow = TRUE), "M")
  d3 <- sc_fixture(matrix(rep(c(0, 1, 2, 1, 0), 3), nrow = 3, byrow = TRUE), "F")
  expect_error(gl.sim.crosses(d3, m5, broodsize = 2, verbose = 0,
                              compliance.check = FALSE),
               "same number of individuals")
})

test_that("mismatched locus panels are refused", {
  # [approved F7] was a DEFECT PIN: the father's loci were ignored and the
  # mother's locus names silently adopted.
  mums <- sc_fixture(matrix(c(0, 1, 2, 1), nrow = 1), "M")
  dads <- sc_fixture(matrix(c(0, 1, 2, 1), nrow = 1), "F")
  locNames(dads) <- c("Z1", "Z2", "Z3", "Z4")
  expect_error(gl.sim.crosses(dads, mums, broodsize = 2, verbose = 0,
                              compliance.check = FALSE),
               "not the same loci")
  # error.check = FALSE opts out of the panel comparison only
  out <- sc_run(dads, mums, broodsize = 2, verbose = 0, error.check = FALSE,
                compliance.check = FALSE)
  expect_equal(locNames(out), c("L1", "L2", "L3", "L4"))
  # unequal locus counts are fatal regardless of error.check
  short <- sc_fixture(matrix(c(0, 1), nrow = 1), "F")
  expect_error(gl.sim.crosses(short, mums, broodsize = 2, verbose = 0,
                              error.check = FALSE, compliance.check = FALSE),
               "same loci")
})

# ---- returned object -------------------------------------------------------

test_that("the returned object is a compliant dartR genlight when compliance.check = TRUE", {
  out1 <- capture.output(males <- gl.keep.pop(testset.gl, pop.list = "Male",
                                             as.pop = "sex", verbose = 0))
  out2 <- capture.output(fems <- gl.keep.pop(testset.gl, pop.list = "Female",
                                             as.pop = "sex", verbose = 0))
  set.seed(32)
  fa <- suppressWarnings({ o <- capture.output(
    v <- gl.keep.ind(males, sample(indNames(males), 10), verbose = 0)); v })
  mo <- suppressWarnings({ o <- capture.output(
    v <- gl.keep.ind(fems, sample(indNames(fems), 10), verbose = 0)); v })
  set.seed(33)
  off <- sc_run(fa, mo, broodsize = 3, n = 10, verbose = 0)
  expect_equal(nInd(off), 10)               # [approved F2] n is honoured
  expect_equal(nLoc(off), nLoc(mo))
  expect_true(all(ploidy(off) == 2))
  expect_equal(locNames(off), locNames(mo))
  expect_equal(indNames(off), paste0("Po_", 1:10))
  expect_equal(nrow(off@other$ind.metrics), nInd(off))
  expect_equal(nrow(off@other$loc.metrics), nLoc(off))
  expect_true(inherits(off@other$ind.metrics, "data.frame"))
  expect_equal(nPop(off), 1)
  # [approved F11] parentage is recorded alongside id and sex
  expect_equal(names(off@other$ind.metrics),
               c("id", "sex", "mother", "father"))
  expect_true(all(off@other$ind.metrics$mother %in% indNames(mo)))
  expect_true(all(off@other$ind.metrics$father %in% indNames(fa)))
  # [approved F13] parental sequence-level locus metrics are carried over
  expect_true("TrimmedSequence" %in% names(off@other$loc.metrics))
  expect_false(any(grepl("array", names(off@other$loc.metrics), fixed = TRUE)))
  # [approved F10] history is this call alone, not the internal helper calls
  expect_length(off@other$history, 1)
  expect_equal(as.character(off@other$history[[1]][[1]]), "gl.sim.crosses")
})

test_that("compliance.check = FALSE returns an object dartR can use", {
  # [approved F5] was a DEFECT PIN: `gl2@other$ind.metrics$sex <- sr` on a
  # fresh genlight assigned into NULL, so ind.metrics became a bare list and
  # loc.metrics, the flags and pop were all absent.
  mums <- sc_fixture(matrix(c(2, 1, 2, 1), nrow = 2, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(2, 1, 0, 1), nrow = 2, byrow = TRUE), "F")
  out <- sc_run(dads, mums, broodsize = 2, verbose = 0, compliance.check = FALSE)
  expect_true(inherits(out@other$ind.metrics, "data.frame"))
  expect_equal(nrow(out@other$ind.metrics), nInd(out))
  expect_equal(nrow(out@other$loc.metrics), nLoc(out))
  expect_false(is.null(out@other$loc.metrics.flags))
  expect_equal(as.character(unique(pop(out))), "pop1")
  filtered <- suppressWarnings(gl.filter.callrate(out, threshold = 0.5,
                                                  verbose = 0))
  expect_s4_class(filtered, "dartR")
  expect_equal(nInd(filtered), nInd(out))
})

test_that("the parent objects are returned untouched", {
  mums <- sc_fixture(matrix(c(2, 1, 2, 1), nrow = 2, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(2, 1, 0, 1), nrow = 2, byrow = TRUE), "F")
  m0 <- mums; d0 <- dads
  set.seed(29)
  invisible(sc_run(dads, mums, broodsize = 3, verbose = 0,
                   compliance.check = FALSE))
  expect_identical(mums, m0)
  expect_identical(dads, d0)
})

test_that("output is reproducible under a fixed seed", {
  mums <- sc_fixture(matrix(c(2, 1, 2, 1), nrow = 2, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(2, 1, 0, 1), nrow = 2, byrow = TRUE), "F")
  set.seed(555)
  a <- sc_run(dads, mums, broodsize = 20, verbose = 0, compliance.check = FALSE)
  set.seed(555)
  b <- sc_run(dads, mums, broodsize = 20, verbose = 0, compliance.check = FALSE)
  expect_identical(as.matrix(a), as.matrix(b))
  expect_identical(a@other$ind.metrics$sex, b@other$ind.metrics$sex)
})

# ---- datatype dispatch ------------------------------------------------------

test_that("SilicoDArT input is refused", {
  # [approved F4] was a DEFECT PIN: utils.check.datatype was called without
  # accept = "SNP", so presence/absence data was crossed as though 1 meant
  # heterozygote and came back scored 0/1/2 at ploidy 2.
  fa <- suppressWarnings({ o <- capture.output(
    v <- gl.keep.ind(testset.gs, indNames(testset.gs)[1:5], verbose = 0)); v })
  mo <- suppressWarnings({ o <- capture.output(
    v <- gl.keep.ind(testset.gs, indNames(testset.gs)[6:10], verbose = 0)); v })
  expect_true(all(ploidy(fa) == 1))
  expect_error(gl.sim.crosses(fa, mo, broodsize = 2, verbose = 0),
               "SilicoDArT")
})

test_that("FBM-backed parents are accepted (by densification)", {
  skip_if_not(exists("gl.gen2fbm"), "gl.gen2fbm not available")
  out1 <- capture.output(males <- gl.keep.pop(testset.gl, pop.list = "Male",
                                             as.pop = "sex", verbose = 0))
  set.seed(43)
  fa <- suppressWarnings({ o <- capture.output(
    v <- gl.keep.ind(males, sample(indNames(males), 6), verbose = 0)); v })
  mo <- suppressWarnings({ o <- capture.output(
    v <- gl.keep.ind(males, sample(indNames(males), 6), verbose = 0)); v })
  ffa <- suppressWarnings({ o <- capture.output(v <- gl.gen2fbm(fa, verbose = 0)); v })
  fmo <- suppressWarnings({ o <- capture.output(v <- gl.gen2fbm(mo, verbose = 0)); v })
  off <- sc_run(ffa, fmo, broodsize = 2, verbose = 0)
  expect_equal(nInd(off), 12)
  expect_true(all(ploidy(off) == 2))
})

# ---- verbosity --------------------------------------------------------------

test_that("verbose = 0 is silent", {
  # [approved F9] was a DEFECT PIN: `cat(report("father --"))`,
  # `cat(report("mother --"))` and the brood-total warning all printed at
  # verbose = 0.
  mums <- sc_fixture(matrix(c(2, 1, 2, 1), nrow = 2, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(2, 1, 0, 1), nrow = 2, byrow = TRUE), "F")
  txt <- capture.output(suppressWarnings(
    res <- gl.sim.crosses(dads, mums, broodsize = 2, n = 1, verbose = 0,
                   compliance.check = FALSE)))
  expect_length(txt, 0)
  txt2 <- capture.output(suppressWarnings(
    res <- gl.sim.crosses(dads, mums, broodsize = 2, verbose = 0,
                   error.check = FALSE, compliance.check = TRUE)))
  expect_length(txt2, 0)
})

test_that("verbose >= 1 reports start and completion", {
  mums <- sc_fixture(matrix(c(2, 1, 2, 1), nrow = 2, byrow = TRUE), "M")
  dads <- sc_fixture(matrix(c(2, 1, 0, 1), nrow = 2, byrow = TRUE), "F")
  txt <- capture.output(suppressWarnings(
    res <- gl.sim.crosses(dads, mums, broodsize = 2, verbose = 2,
                   compliance.check = FALSE)))
  expect_true(any(grepl("Starting gl.sim.crosses", txt)))
  expect_true(any(grepl("Completed: gl.sim.crosses", txt)))
  # [approved F9] the preamble always runs, even with error.check = FALSE
  txt2 <- capture.output(suppressWarnings(
    res <- gl.sim.crosses(dads, mums, broodsize = 2, verbose = 2,
                   error.check = FALSE, compliance.check = FALSE)))
  expect_true(any(grepl("Starting gl.sim.crosses", txt2)))
  expect_true(any(grepl("Completed: gl.sim.crosses", txt2)))
})
