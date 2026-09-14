# Characterization baseline for gl.compliance.check (function-review campaign,
# Phase A). Pins CURRENT behaviour, including known defects (marked BUG); a
# behaviour change that breaks a BUG pin should map to an approved finding in
# function-review/reports/dartR.base/gl.compliance.check.md.
#
# Position-slot expectations follow PR #330 (position-genome-only), which is
# merged into the local integration branch: @position is never filled from
# loc.metrics$SnpPosition, and a stale copy is cleared.

# ---- helpers ---------------------------------------------------------------

make_gl <- function(m, ploidy = 2, pop = "a") {
  g <- new("genlight", gen = m, ploidy = ploidy)
  indNames(g) <- paste0("i", seq_len(nrow(m)))
  locNames(g) <- paste0("L", seq_len(ncol(m)))
  pop(g) <- factor(rep(pop, nrow(m)))
  g
}

# ---- baseline: happy path --------------------------------------------------

test_that("baseline run on testset2.gl is silent at verbose=0 and idempotent", {
  out <- capture.output(a <- gl.compliance.check(testset2.gl, verbose = 0))
  expect_length(out, 0)
  expect_s4_class(a, "dartR")
  expect_equal(nLoc(a), nLoc(testset2.gl))
  expect_equal(nInd(a), nInd(testset2.gl))
  # flags exist and the calculable metrics are flagged TRUE after recalc
  fl <- a@other$loc.metrics.flags
  expect_true(all(unlist(fl[c("AvgPIC", "CallRate", "maf", "FreqHets",
                              "FreqHomRef", "FreqHomSnp")])))
  # testset2.gl contains monomorphic and all-NA loci; flags record that
  expect_false(fl$monomorphs)
  expect_false(fl$allna)
  # second run changes nothing but history
  out2 <- capture.output(b <- gl.compliance.check(a, verbose = 0))
  expect_length(out2, 0)
  a@other$history <- b@other$history <- list()
  expect_true(isTRUE(all.equal(a, b)))
})

test_that("BUG(F9): one call appends two history entries (internal leak)", {
  a <- gl.compliance.check(testset2.gl, verbose = 0)
  added <- length(a@other$history) - length(testset2.gl@other$history)
  expect_equal(added, 2)  # gl.recalc.metrics internal call + own match.call
  calls <- vapply(a@other$history, function(h) as.character(h[[1]]), "")
  expect_true("gl.recalc.metrics" %in% calls)
  expect_true("gl.compliance.check" %in% calls)
})

test_that("baseline run on testset2.gs (SilicoDArT) is silent at verbose=0", {
  out <- capture.output(a <- gl.compliance.check(testset2.gs, verbose = 0))
  expect_length(out, 0)
  expect_s4_class(a, "dartR")
  expect_true(a@other$loc.metrics.flags$OneRatio)
  expect_true(a@other$loc.metrics.flags$PIC)
})

# ---- repairs that work -----------------------------------------------------

test_that("bare adegenet genlight is repaired: class, pop, metrics, names", {
  set.seed(1)
  m <- matrix(sample(c(0L, 1L, 2L, NA), 40, TRUE), nrow = 4)
  g <- new("genlight", gen = m, ploidy = 2)  # no @other, no pop, no names
  out <- capture.output(res <- gl.compliance.check(g, verbose = 0))
  expect_length(out, 0)
  expect_s4_class(res, "dartR")
  expect_equal(levels(pop(res)), "pop1")
  expect_equal(locNames(res), paste0("Loc", 1:10))
  expect_true(is.data.frame(res@other$loc.metrics))
  expect_equal(nrow(res@other$loc.metrics), 10)
  expect_equal(res@other$ind.metrics$id, indNames(res))
  expect_equal(res@other$loc.metrics$CallRate,
               unname(colMeans(!is.na(as.matrix(res)))))
})

test_that("NULL ind.metrics is repaired with an id column", {
  g <- testset2.gl[1:10, 1:20]
  g@other$ind.metrics <- NULL
  res <- gl.compliance.check(g, verbose = 0)
  expect_equal(names(res@other$ind.metrics), "id")
  expect_equal(res@other$ind.metrics$id, indNames(res))
})

test_that("empty loc.all is filled with the dummy A/C", {
  g <- testset2.gl[1:10, 1:20]
  g@loc.all <- NULL
  res <- gl.compliance.check(g, verbose = 0)
  expect_equal(unique(res@loc.all), "A/C")
})

test_that("latlong/long misspellings are renamed to latlon/lon", {
  g <- testset2.gl[1:10, 1:20]
  g@other$latlong <- g@other$latlon
  g@other$latlon <- NULL
  res <- gl.compliance.check(g, verbose = 0)
  expect_null(res@other$latlong)
  expect_false(is.null(res@other$latlon))
})

test_that("NULL @ploidy with uniform diploid data is repaired to 2", {
  set.seed(2)
  m <- matrix(sample(c(0L, 1L, 2L), 40, TRUE), nrow = 4)
  m[1, 1] <- 2L  # every individual carries a 2 so computed ploidy is uniform
  m[2, 1] <- 2L; m[3, 2] <- 2L; m[4, 2] <- 2L
  g <- make_gl(m)
  g@ploidy <- NULL
  res <- gl.compliance.check(g, verbose = 0)
  expect_true(all(ploidy(res) == 2))
})

test_that("duplicated individual names are made unique", {
  g <- testset2.gl[1:10, 1:20]
  indNames(g) <- c("dup", "dup", paste0("i", 3:10))
  g@other$ind.metrics$id <- indNames(g)
  res <- gl.compliance.check(g, verbose = 0)
  expect_false(any(duplicated(indNames(res))))
  expect_equal(indNames(res)[2], "dup_1")
})

# ---- position slot (post-#330 convention) ----------------------------------

test_that("@position handling is consistent with this tree's position convention", {
  # State-agnostic: pre-#330 code fills @position from loc.metrics$SnpPosition;
  # post-#330 (position-genome-only) leaves it NULL for genome coordinates.
  g <- testset2.gl[1:10, 1:20]
  g@position <- NULL
  res <- gl.compliance.check(g, verbose = 0)
  if (is.null(res@position)) {
    # post-#330 convention
    expect_null(res@position)  # stays NULL: reserved for genome coordinates
    g2 <- testset2.gl[1:10, 1:20]
    g2@position <- as.integer(g2@other$loc.metrics$SnpPosition)
    res2 <- gl.compliance.check(g2, verbose = 0)
    expect_null(res2@position)  # provable stale copy is cleared
  } else {
    # pre-#330 behaviour: filled from the tag offset
    expect_equal(as.integer(res@position),
                 as.integer(g@other$loc.metrics$SnpPosition))
  }
})

# ---- known defects, pinned as-is (BUG) -------------------------------------

test_that("F1 fixed: out-of-range genotype is reported even when NAs present", {
  # [approved F1 minimal] was BUG(F1): max(mat) is NA and NA %in%
  # c(0,1,2,NA) is TRUE, so a genotype of 5 printed "confirmed". The check
  # now tests the value set exactly; gating and warn-not-stop are retained
  # by decision (2026-09-05), so the 5 is reported, not repaired.
  m <- matrix(c(0L, 1L, 2L, 5L, 0L, NA, 1L, 2L, 0L, 1L, 2L, 0L), nrow = 4)
  g <- make_gl(m)
  out <- capture.output(res <- gl.compliance.check(g, verbose = 2))
  expect_false(any(grepl("scored NA, 0, 1 or 2 confirmed", out)))
  expect_true(any(grepl("must be scored NA, 0 or 1 or 2", out)))
  expect_s4_class(res, "dartR")                       # still warn-not-stop
  expect_equal(max(as.matrix(res), na.rm = TRUE), 5)  # still no repair
})

test_that("F1 retained semantics: coding error is printed but not raised, silent at v0", {
  # [approved F1 minimal] warn-not-stop and verbose >= 1 gating kept as-is;
  # only the NA-vacuity was fixed. This path behaved the same before the fix.
  m <- matrix(c(0L, 1L, 2L, 5L, 0L, 1L, 1L, 2L, 0L, 1L, 2L, 0L), nrow = 4)
  g <- make_gl(m)
  out <- capture.output(res <- gl.compliance.check(g, verbose = 2))
  expect_true(any(grepl("must be scored NA, 0 or 1 or 2", out)))
  expect_s4_class(res, "dartR")   # returns anyway, no stop()
  out0 <- capture.output(res0 <- gl.compliance.check(g, verbose = 0))
  expect_length(out0, 0)          # coding violation invisible at verbose=0
})

test_that("F2 fixed: all-monomorphic input completes with monomorphs flag FALSE", {
  # [approved F2 minimal] was BUG(F2): gl.filter.monomorphs -> gl.drop.loc
  # errored "Subsetting resulted in zero loci". That error is now caught and
  # read as the all-monomorphic case: flag-only, no crash, no repair.
  m <- matrix(c(0L, 0L, 0L, 0L, 2L, 2L, 2L, 2L, 0L, 0L, 0L, NA), nrow = 4)
  g <- make_gl(m)
  out <- capture.output(res <- gl.compliance.check(g, verbose = 0))
  expect_length(out, 0)
  expect_s4_class(res, "dartR")
  expect_equal(nLoc(res), 3)
  expect_false(res@other$loc.metrics.flags$monomorphs)
  out2 <- capture.output(res2 <- gl.compliance.check(g, verbose = 2))
  expect_true(any(grepl("All loci are monomorphic", out2)))
})

test_that("F2 fixed: all-all-NA input completes at verbose 0/1 with flags FALSE", {
  # [approved F2 minimal] the same filter-to-count idiom crashed via
  # gl.filter.allna when every locus was all missing; now caught the same
  # way. At verbose >= 2 the SAME root cause (gl.filter.allna cannot
  # return zero loci) still crashes earlier, inside utils.check.datatype's
  # verbose>1 courtesy check -- pre-existing at ddaed27, outside this
  # function's file, recorded in the review report.
  m <- matrix(NA_integer_, nrow = 4, ncol = 3)
  g <- make_gl(m)
  out <- capture.output(res <- gl.compliance.check(g, verbose = 0))
  expect_length(out, 0)
  expect_s4_class(res, "dartR")
  expect_false(res@other$loc.metrics.flags$allna)
  expect_false(res@other$loc.metrics.flags$monomorphs)
})

test_that("BUG(F3): NULL @ploidy with heterogeneous computed ploidy crashes", {
  m <- rbind(c(0L, 1L, 1L, 0L),   # per-individual max 1 -> computed ploidy 1
             c(0L, 1L, 0L, 1L),
             c(2L, 1L, 0L, 2L),   # computed ploidy 2
             c(0L, 2L, 1L, 0L))
  g <- make_gl(m)
  g@ploidy <- NULL
  expect_error(gl.compliance.check(g, verbose = 0),
               "condition has length > 1")
})

test_that("BUG(F4): make.unique repair desynchronises ind.metrics$id", {
  g <- testset2.gl[1:10, 1:20]
  indNames(g) <- c("dup", "dup", paste0("i", 3:10))
  g@other$ind.metrics$id <- indNames(g)
  res <- gl.compliance.check(g, verbose = 0)
  expect_false(identical(indNames(res), as.character(res@other$ind.metrics$id)))
})

test_that("BUG(F5): loc.metrics row desync crashes opaquely, designed warning unreachable", {
  g <- testset2.gl[1:10, 1:20]
  g@other$loc.metrics <- g@other$loc.metrics[1:10, ]
  expect_error(gl.compliance.check(g, verbose = 0),
               "replacement has")
})

test_that("BUG(F6): ind.metrics row desync passes silently", {
  g <- testset2.gl[1:10, 1:20]
  g@other$ind.metrics <- g@other$ind.metrics[1:5, ]
  out <- capture.output(res <- gl.compliance.check(g, verbose = 0))
  expect_length(out, 0)
  expect_false(nrow(res@other$ind.metrics) == nInd(res))
})

test_that("BUG(F7): an all-NA individual is retained and allna flag reads TRUE", {
  set.seed(7)
  m <- matrix(sample(c(0L, 1L, 2L), 60, TRUE), nrow = 5)
  m[1, ] <- NA
  g <- make_gl(m)
  res <- gl.compliance.check(g, verbose = 0)
  expect_equal(nInd(res), 5)
  expect_true(res@other$loc.metrics.flags$allna)
})

test_that("BUG(F8): non-uniform stamped ploidy passes without detection", {
  g <- testset2.gl[1:10, 1:20]
  g@ploidy <- as.integer(c(3, 3, rep(2, 8)))
  out <- capture.output(res <- gl.compliance.check(g, verbose = 0))
  expect_length(out, 0)
  expect_equal(sort(unique(ploidy(res))), c(2L, 3L))
})

test_that("BUG(F10): an NA locus name crashes in metric recalculation", {
  g <- testset2.gl[1:10, 1:20]
  ln <- locNames(g)
  ln[5] <- NA
  locNames(g) <- ln
  expect_error(gl.compliance.check(g, verbose = 0),
               "missing values")
})

test_that("BUG(F11): duplicated locus names pass silently", {
  g <- testset2.gl[1:10, 1:20]
  locNames(g) <- c("dupL", "dupL", paste0("L", 3:20))
  res <- gl.compliance.check(g, verbose = 0)
  expect_true(any(duplicated(locNames(res))))
})
