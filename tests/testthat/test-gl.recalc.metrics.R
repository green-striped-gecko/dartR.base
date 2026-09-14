# Characterization tests for gl.recalc.metrics
# Baseline snapshotted before review (dev at ddaed27; helper battery at
# integration-local ed99203). These assertions pin behaviour against the
# review in function-review/reports/dartR.base/gl.recalc.metrics.md.
#
# 2026-09-08: findings F2 to F10 were approved and applied. The pins they
# encoded are flipped and annotated [approved Fn] with the report's finding
# IDs. The pin tags written in the first draft of this file used a numbering
# that differs from the report for four items -- file F2 = report F3
# (history), file F3 = report F2 (loc.metrics), file F6 = report F9 (mono.rm
# validation), file F7 = report F6 (verbose) -- so the annotations below
# carry the report ID.
#
# The [pins defect F1] assertions (stale rdepth, stale AvgReadDepth) are
# unchanged. F1 is deferred (custodian, 2026-09-08); they remain the record
# of that defect and must keep passing.

sub_snp <- function() {
  gl.drop.pop(testset.gl, pop.list = popNames(testset.gl)[1:10], verbose = 0)
}

sub_silico <- function() {
  gl.drop.pop(testset.gs, pop.list = popNames(testset.gs)[1:5], verbose = 0)
}

# ---------------------------------------------------------------------------
# 1. Metric correctness (SNP)
# ---------------------------------------------------------------------------

test_that("every SNP metric matches an independent recomputation from as.matrix", {
  x <- sub_snp()
  r <- gl.recalc.metrics(x, verbose = 0)
  t <- as.matrix(r)
  c0 <- colSums(t == 0, na.rm = TRUE)
  c1 <- colSums(t == 1, na.rm = TRUE)
  c2 <- colSums(t == 2, na.rm = TRUE)
  ctot <- c0 + c1 + c2
  orr <- (c0 + c1) / ctot
  ors <- (c1 + c2) / ctot
  picr <- 1 - (orr^2 + (1 - orr)^2)
  pics <- 1 - (ors^2 + (1 - ors)^2)
  alf <- colMeans(t, na.rm = TRUE) / 2
  lm <- r@other$loc.metrics

  expect_equal(unname(lm$CallRate),
               unname(signif(1 - colSums(is.na(t)) / nInd(r), 6)))
  expect_equal(unname(lm$OneRatioRef), unname(orr))
  expect_equal(unname(lm$OneRatioSnp), unname(ors))
  expect_equal(unname(lm$PICRef), unname(picr))
  expect_equal(unname(lm$PICSnp), unname(pics))
  expect_equal(unname(lm$AvgPIC), unname((picr + pics) / 2))
  expect_equal(unname(lm$FreqHomRef), unname(c0 / ctot))
  expect_equal(unname(lm$FreqHomSnp), unname(c2 / ctot))
  expect_equal(unname(lm$FreqHets), unname(c1 / ctot))
  expect_equal(unname(lm$maf), unname(pmin(alf, 1 - alf)))
})

test_that("platypus.gl gives the same exact agreement (1000 loci)", {
  r <- gl.recalc.metrics(platypus.gl, verbose = 0)
  t <- as.matrix(r)
  c0 <- colSums(t == 0, na.rm = TRUE)
  c1 <- colSums(t == 1, na.rm = TRUE)
  c2 <- colSums(t == 2, na.rm = TRUE)
  ctot <- c0 + c1 + c2
  expect_equal(unname(r@other$loc.metrics$FreqHets), unname(c1 / ctot))
  expect_equal(unname(r@other$loc.metrics$CallRate),
               unname(signif(1 - colSums(is.na(t)) / nInd(r), 6)))
})

# ---------------------------------------------------------------------------
# 2. Staleness matrix: which loc.metrics columns are refreshed
# ---------------------------------------------------------------------------

test_that("the refreshed column set is exactly the ten SNP metrics", {
  x <- sub_snp()
  r <- gl.recalc.metrics(x, verbose = 0)
  pre <- x@other$loc.metrics
  post <- r@other$loc.metrics
  changed <- names(post)[vapply(names(post), function(cn) {
    !is.null(pre[[cn]]) && !isTRUE(all.equal(pre[[cn]], post[[cn]]))
  }, logical(1))]
  expect_setequal(changed,
                  c("CallRate", "OneRatioRef", "OneRatioSnp", "PICRef",
                    "PICSnp", "AvgPIC", "FreqHomRef", "FreqHomSnp",
                    "FreqHets", "maf"))
})

test_that("[pins defect F1] rdepth is left stale after individuals are dropped", {
  x <- sub_snp()
  r <- gl.recalc.metrics(x, verbose = 0)
  # rdepth is unchanged by the recalculation ...
  expect_equal(x@other$loc.metrics$rdepth, r@other$loc.metrics$rdepth)
  # ... although its defining formula (gl.read.dart) uses the very metrics
  # that were just refreshed.
  lm <- r@other$loc.metrics
  implied <- round(lm$OneRatioRef * lm$AvgCountRef +
                     lm$OneRatioSnp * lm$AvgCountSnp, 1)
  ok <- !is.na(implied)
  # The formula is exact on the unsubsetted dataset:
  lm0 <- testset.gl@other$loc.metrics
  expect_equal(as.numeric(lm0$rdepth),
               as.numeric(round(lm0$OneRatioRef * lm0$AvgCountRef +
                                  lm0$OneRatioSnp * lm0$AvgCountSnp, 1)))
  # After subsetting, the stored value disagrees for most loci:
  expect_gt(sum(abs(lm$rdepth[ok] - implied[ok]) > 0.05), 100)
})

test_that("[pins defect F1] SilicoDArT AvgReadDepth is left stale", {
  s <- sub_silico()
  r <- gl.recalc.metrics(s, verbose = 0)
  expect_equal(s@other$loc.metrics$AvgReadDepth,
               r@other$loc.metrics$AvgReadDepth)
})

# ---------------------------------------------------------------------------
# 3. Flags (DAT4)
# ---------------------------------------------------------------------------

test_that("SNP call sets the ten metric flags TRUE", {
  x <- sub_snp()
  r <- gl.recalc.metrics(x, verbose = 0)
  fl <- r@other$loc.metrics.flags
  for (f in c("AvgPIC", "OneRatioRef", "OneRatioSnp", "PICRef", "PICSnp",
              "CallRate", "maf", "FreqHets", "FreqHomRef", "FreqHomSnp")) {
    expect_true(fl[[f]], info = f)
  }
  expect_false(fl$monomorphs)
})

test_that("[approved F4] a stale monomorphs = TRUE claim is re-examined", {
  x <- sub_snp()
  x@other$loc.metrics.flags$monomorphs <- TRUE
  r <- gl.recalc.metrics(x, verbose = 0)
  # the object still holds monomorphic and all-NA loci, so the flag must not
  # claim otherwise
  expect_false(r@other$loc.metrics.flags$monomorphs)
  # and an object that genuinely has none gets TRUE from the same check
  clean <- gl.filter.monomorphs(x, verbose = 0)
  expect_true(gl.recalc.metrics(clean,
                                verbose = 0)@other$loc.metrics.flags$monomorphs)
  # the check agrees with the gl.filter.monomorphs definition of monomorphic
  m <- as.matrix(x)
  n_mono <- sum(apply(m, 2, function(g) {
    all(g == 0, na.rm = TRUE) || all(g == 2, na.rm = TRUE) || all(is.na(g))
  }))
  expect_gt(n_mono, 0)
})

test_that("SilicoDArT sets only CallRate, OneRatio and PIC", {
  s <- sub_silico()
  r <- gl.recalc.metrics(s, verbose = 0)
  fl <- r@other$loc.metrics.flags
  expect_true(fl$CallRate)
  expect_true(fl$OneRatio)
  expect_true(fl$PIC)
  expect_false(fl$maf)
  expect_false(fl$FreqHets)
  expect_false(fl$AvgPIC)
})

# ---------------------------------------------------------------------------
# 4. Datatype dispatch
# ---------------------------------------------------------------------------

test_that("SilicoDArT metrics match an independent recomputation", {
  s <- sub_silico()
  r <- gl.recalc.metrics(s, verbose = 0)
  t <- as.matrix(r)
  onerat <- colMeans(t == 1, na.rm = TRUE)
  expect_equal(unname(r@other$loc.metrics$OneRatio), unname(onerat))
  expect_equal(unname(r@other$loc.metrics$PIC),
               unname(1 - (onerat^2 + (1 - onerat)^2)))
  expect_equal(unname(r@other$loc.metrics$CallRate),
               unname(signif(1 - colSums(is.na(t)) / nInd(r), 6)))
  # no SNP-only columns are invented on presence/absence data
  expect_setequal(names(r@other$loc.metrics), names(s@other$loc.metrics))
})

# ---------------------------------------------------------------------------
# 5. History (FS8) and the return contract
# ---------------------------------------------------------------------------

test_that("a direct call appends exactly one history entry recording the call", {
  x <- sub_snp()
  n <- length(x@other$history)
  r <- gl.recalc.metrics(x, verbose = 0)
  expect_equal(length(r@other$history), n + 1L)
  expect_equal(as.character(r@other$history[[n + 1L]][[1]]),
               "gl.recalc.metrics")
})

test_that("[approved F3] mono.rm = TRUE appends one entry, not two", {
  x <- sub_snp()
  n <- length(x@other$history)
  r <- gl.recalc.metrics(x, mono.rm = TRUE, verbose = 0)
  expect_equal(length(r@other$history), n + 1L)
  expect_equal(as.character(r@other$history[[n + 1L]][[1]]),
               "gl.recalc.metrics")
  # the internal gl.filter.monomorphs step is not recorded
  expect_false("gl.filter.monomorphs" %in%
                 vapply(r@other$history, function(e) as.character(e[[1]]), ""))
})

test_that("[approved F3] a nested call appends nothing; compliance checks stop multiplying", {
  n <- length(testset.gl@other$history)
  b <- gl.compliance.check(testset.gl, verbose = 0)
  expect_equal(length(b@other$history), n + 1L)
  expect_equal(as.character(b@other$history[[n + 1L]][[1]]),
               "gl.compliance.check")
  b2 <- gl.compliance.check(b, verbose = 0)
  expect_equal(length(b2@other$history), n + 2L)
  b3 <- gl.compliance.check(b2, verbose = 0)
  expect_equal(length(b3@other$history), n + 3L)
  # no gl.recalc.metrics entry anywhere in the chain
  expect_false("gl.recalc.metrics" %in%
                 vapply(b3@other$history, function(e) as.character(e[[1]]), ""))
})

test_that("the return is by value: discarding it leaves the caller's object untouched", {
  x <- sub_snp()
  before_cr <- x@other$loc.metrics$CallRate
  before_h <- length(x@other$history)
  invisible(gl.recalc.metrics(x, verbose = 0))
  expect_identical(before_cr, x@other$loc.metrics$CallRate)
  expect_equal(length(x@other$history), before_h)
})

test_that("repeated calls are idempotent in metrics but not in history", {
  x <- sub_snp()
  a <- gl.recalc.metrics(x, verbose = 0)
  b <- gl.recalc.metrics(a, verbose = 0)
  expect_equal(a@other$loc.metrics, b@other$loc.metrics)
  expect_equal(length(b@other$history), length(a@other$history) + 1L)
})

# ---------------------------------------------------------------------------
# 6. Genotypes and object identity
# ---------------------------------------------------------------------------

test_that("genotypes, individuals, populations and ploidy are untouched", {
  x <- sub_snp()
  r <- gl.recalc.metrics(x, verbose = 0)
  expect_identical(as.matrix(x), as.matrix(r))
  expect_identical(indNames(x), indNames(r))
  expect_identical(pop(x), pop(r))
  expect_identical(ploidy(x), ploidy(r))
  expect_identical(x@other$ind.metrics, r@other$ind.metrics)
})

# ---------------------------------------------------------------------------
# 7. mono.rm
# ---------------------------------------------------------------------------

test_that("mono.rm = TRUE drops monomorphic and all-NA loci and keeps metadata in step", {
  x <- sub_snp()
  allna_before <- sum(colSums(!is.na(as.matrix(x))) == 0)
  expect_gt(allna_before, 0)
  r <- gl.recalc.metrics(x, mono.rm = TRUE, verbose = 0)
  expect_lt(nLoc(r), nLoc(x))
  expect_equal(nrow(r@other$loc.metrics), nLoc(r))   # DAT2
  expect_equal(sum(colSums(!is.na(as.matrix(r))) == 0), 0)
  expect_true(r@other$loc.metrics.flags$monomorphs)
})

test_that("[approved F5] mono.rm = TRUE on all-monomorphic data completes with a message", {
  x <- sub_snp()
  m <- as.matrix(x)
  const <- apply(m, 2, function(g) length(unique(stats::na.omit(g))) <= 1)
  mm <- x[, which(const)[1:5]]
  expect_silent(r <- gl.recalc.metrics(mm, mono.rm = TRUE, verbose = 0))
  expect_equal(nLoc(r), 5)
  expect_false(r@other$loc.metrics.flags$monomorphs)
  o <- capture.output(r <- gl.recalc.metrics(mm, mono.rm = TRUE, verbose = 1))
  expect_true(any(grepl("monomorphic", o)))
  expect_false(any(grepl("Subsetting resulted in zero loci", o)))
})

test_that("[approved F9] mono.rm is validated with a dartR error", {
  x <- sub_snp()
  expect_error(gl.recalc.metrics(x, mono.rm = "yes", verbose = 0),
               "mono.rm must be a single logical value")
  expect_error(gl.recalc.metrics(x, mono.rm = NA, verbose = 0),
               "mono.rm must be a single logical value")
  expect_error(gl.recalc.metrics(x, mono.rm = NULL, verbose = 0),
               "mono.rm must be a single logical value")
  expect_error(gl.recalc.metrics(x, mono.rm = c(TRUE, TRUE), verbose = 0),
               "mono.rm must be a single logical value")
})

# ---------------------------------------------------------------------------
# 8. Edge cases and object-shape guards
# ---------------------------------------------------------------------------

test_that("[approved F2] a missing loc.metrics slot is reported and rebuilt, never fabricated", {
  x <- sub_snp()

  # Many loci: a conforming table is built, no "replacement has N rows" error.
  many <- x
  many@other$loc.metrics <- NULL
  expect_silent(r <- gl.recalc.metrics(many, verbose = 0))
  expect_true(is.data.frame(r@other$loc.metrics))
  expect_equal(nrow(r@other$loc.metrics), nLoc(r))
  o <- capture.output(r <- gl.recalc.metrics(many, verbose = 1))
  expect_true(any(grepl("no locus metrics data frame found", o)))

  # One locus: the flags table is no longer reached by `$` partial matching,
  # so it is not returned as the locus metrics table.
  one <- x[, 1]
  one@other$loc.metrics <- NULL
  r <- gl.recalc.metrics(one, verbose = 0)
  expect_true(is.data.frame(r@other$loc.metrics))
  expect_equal(nrow(r@other$loc.metrics), 1)
  expect_false(identical(r@other$loc.metrics, r@other$loc.metrics.flags))
  expect_false("monomorphs" %in% names(r@other$loc.metrics))
  expect_true(all(c("CallRate", "AvgPIC") %in% names(r@other$loc.metrics)))

  # A table whose rows do not track the loci is not repairable here.
  bad <- x
  bad@other$loc.metrics <- bad@other$loc.metrics[1:10, , drop = FALSE]
  expect_error(gl.recalc.metrics(bad, verbose = 0),
               "must track loci one for one")
})

test_that("[approved F2] a plain genlight gets loc.metrics as a data frame", {
  x <- sub_snp()
  gg <- new("genlight", gen = as.matrix(x), ploidy = 2)
  r <- gl.recalc.metrics(gg, verbose = 0)
  expect_true(is.data.frame(r@other$loc.metrics))   # DAT2: rows track loci
  expect_equal(nrow(r@other$loc.metrics), nLoc(r))
  # and the table is subsettable row-wise, as every gl.filter.* needs
  expect_equal(nrow(r@other$loc.metrics[1:3, , drop = FALSE]), 3)
})

test_that("missing individual metric columns are recreated", {
  x <- sub_snp()
  x@other$loc.metrics$CallRate <- NULL
  x@other$loc.metrics$AvgPIC <- NULL
  r <- gl.recalc.metrics(x, verbose = 0)
  expect_false(is.null(r@other$loc.metrics$CallRate))
  expect_false(is.null(r@other$loc.metrics$AvgPIC))
})

test_that("a missing flags slot is tolerated and rebuilt", {
  x <- sub_snp()
  x@other$loc.metrics.flags <- NULL
  r <- gl.recalc.metrics(x, verbose = 0)
  expect_true(r@other$loc.metrics.flags$CallRate)
})

test_that("all-NA loci yield NaN metrics without a warning at verbose 0", {
  x <- sub_snp()
  r <- gl.recalc.metrics(x, verbose = 0)
  i <- which(colSums(!is.na(as.matrix(r))) == 0)[1]
  expect_false(is.na(i))
  expect_equal(r@other$loc.metrics$CallRate[i], 0)
  expect_true(is.nan(r@other$loc.metrics$AvgPIC[i]))
  expect_true(is.nan(r@other$loc.metrics$FreqHets[i]))
})

test_that("a single-locus object is handled", {
  x <- sub_snp()
  one <- x[, 1]
  one@other$loc.metrics <- x@other$loc.metrics[1, , drop = FALSE]
  r <- gl.recalc.metrics(one, verbose = 0)
  expect_equal(nLoc(r), 1)
  expect_equal(nrow(r@other$loc.metrics), 1)
})

# ---------------------------------------------------------------------------
# 9. Verbosity
# ---------------------------------------------------------------------------

test_that("verbose = 0 is fully silent on both streams and in both branches", {
  x <- sub_snp()
  expect_equal(length(capture.output(r <- gl.recalc.metrics(x, verbose = 0))), 0)
  expect_equal(length(capture.output(r <- gl.recalc.metrics(x, mono.rm = TRUE,
                                                            verbose = 0))), 0)
  expect_equal(length(capture.output(r <- gl.recalc.metrics(x, verbose = 0),
                                     type = "message")), 0)
  s <- sub_silico()
  expect_equal(length(capture.output(r <- gl.recalc.metrics(s, verbose = 0))), 0)
})

test_that("[approved F6] verbose = 1 prints this function's banners only", {
  x <- sub_snp()
  o <- capture.output(r <- gl.recalc.metrics(x, verbose = 1))
  # VRB1: two lines (begin and end) for this function alone.
  expect_equal(length(o), 2)
  expect_false(any(grepl("utils.recalc", o)))
  # at verbose = 2 the consolidated report replaces the helpers' repeats
  o2 <- capture.output(r <- gl.recalc.metrics(x, verbose = 2))
  expect_false(any(grepl("utils.recalc", o2)))
  expect_equal(sum(grepl("monomorphic", o2)), 1)
  expect_true(any(grepl("Locus metrics recalculated", o2)))
})
