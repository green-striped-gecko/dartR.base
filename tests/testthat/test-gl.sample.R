# Characterization tests for gl.sample
#
# Captured at package commit ddaed27 (upstream/dev) as a baseline, then updated
# in the same commit as the review fix. Assertions that changed are marked
# "# [approved Fn]" and name the finding ID from
# function-review/reports/dartR.base/gl.sample.md. Every other assertion is
# unchanged from the baseline, so a diff outside the marked lines is a
# regression.

# Reference draw: reproduces gl.sample's own draw so a test can predict which
# source individuals a seed selects. Indexes into the pool rather than passing
# the pool to sample(), which is the F1 fix.
ref_draw <- function(x, nsample, replace = TRUE, onepop = FALSE) {
  if (onepop) {
    pools <- list(seq_len(nInd(x)))
  } else {
    pops <- pop(x)
    sizes <- table(pops)
    sizes <- sizes[sizes > 0]
    pools <- lapply(names(sizes), function(p) which(pops == p))
  }
  unlist(lapply(pools, function(idx)
    idx[sample.int(length(idx), nsample, replace = replace)]))
}

# ---------------------------------------------------------------------------
# Sampling correctness: counts, pool, replacement
# ---------------------------------------------------------------------------

test_that("gl.sample draws nsample individuals from every population", {
  x <- platypus.gl                      # 81 ind, 1000 loc, 3 pops, no singletons
  set.seed(1)
  out <- gl.sample(x, nsample = 4, replace = TRUE, verbose = 0)

  expect_s4_class(out, "dartR")
  expect_equal(nInd(out), 4 * nPop(x))
  expect_equal(nLoc(out), nLoc(x))
  expect_true(all(table(as.character(pop(out))) == 4))
  expect_true(all(ploidy(out) == 2))
})

# Output names carry the NN_ prefix gl.sample adds. Strip it to recover the
# source individual.
src_names <- function(x) sub("^[0-9]+_", "", indNames(x))

test_that("gl.sample draws only from the requested population's own members", {
  x <- platypus.gl
  set.seed(2)
  out <- gl.sample(x, nsample = 5, replace = TRUE, verbose = 0)

  src <- src_names(out)
  expect_true(all(src %in% indNames(x)))
  src_pop <- as.character(pop(x))[match(src, indNames(x))]
  expect_equal(src_pop, as.character(pop(out)))
})

test_that("gl.sample replace = FALSE draws each individual at most once", {
  x <- platypus.gl
  set.seed(3)
  out <- gl.sample(x, nsample = 10, replace = FALSE, verbose = 0)

  src <- src_names(out)
  expect_false(any(duplicated(src)))
  expect_equal(nInd(out), 10 * nPop(x))
})

test_that("gl.sample replace = TRUE can draw the same individual more than once", {
  x <- platypus.gl
  # 81 individuals drawn 20 times with replacement: collisions are near certain
  set.seed(4)
  out <- gl.sample(x, nsample = 20, replace = TRUE, onepop = TRUE, verbose = 0)
  src <- src_names(out)
  expect_true(any(duplicated(src)))
  # [approved F2] the zero-padded ordinal prefix is what keeps the returned
  # names unique; rbind's __indN suffix is gone with the rbind assembly
  expect_false(any(duplicated(indNames(out))))
  expect_false(any(grepl("__ind[0-9]+$", indNames(out))))
})

test_that("gl.sample is reproducible under a fixed seed", {
  x <- platypus.gl
  set.seed(11); a <- gl.sample(x, nsample = 4, replace = TRUE, verbose = 0)
  set.seed(11); b <- gl.sample(x, nsample = 4, replace = TRUE, verbose = 0)

  expect_identical(indNames(a), indNames(b))
  expect_identical(as.matrix(a), as.matrix(b))
})

test_that("gl.sample onepop = TRUE ignores population structure", {
  x <- platypus.gl
  set.seed(5)
  out <- gl.sample(x, nsample = 6, onepop = TRUE, replace = TRUE, verbose = 0)

  # nsample is a whole-object total, not a per-population count
  expect_equal(nInd(out), 6)
  # original population labels are restored on the drawn individuals
  expect_true(all(as.character(pop(out)) %in% levels(pop(x))))
})

# ---------------------------------------------------------------------------
# F1: every population is drawn from its own members, singletons included
# ---------------------------------------------------------------------------

test_that("gl.sample draws populations with exactly one member", {
  x <- testset.gl
  singles <- names(which(table(pop(x)) == 1))
  expect_setequal(singles, c("EmmacNormLeic", "EmmacNormSalt"))

  set.seed(42)
  out <- gl.sample(x, nsample = 3, replace = TRUE, verbose = 0)

  # [approved F1] sample(v, n) with length(v) == 1 drew from 1:v, so a
  # single-member population was replaced by unrelated individuals from
  # elsewhere in the object. sample.int(length(v), n) draws from v itself.
  expect_true(all(singles %in% as.character(pop(out))))
  expect_equal(length(unique(as.character(pop(out)))), nPop(x))
  expect_true(all(table(as.character(pop(out))) == 3))
  expect_equal(nInd(out), 3 * nPop(x))
})

test_that("gl.sample returns a singleton population's own member, never a stranger", {
  x <- testset.gl
  singles <- names(which(table(pop(x)) == 1))
  # [approved F1] across 20 seeds each singleton population contributes its own
  # single member, repeated nsample times
  for (s in 1:20) {
    set.seed(s)
    out <- gl.sample(x, nsample = 3, replace = TRUE, verbose = 0)
    src <- src_names(out)
    for (p in singles) {
      own <- indNames(x)[as.character(pop(x)) == p]
      expect_equal(src[as.character(pop(out)) == p], rep(own, 3))
    }
  }
})

test_that("gl.sample assigns every drawn individual its own population label", {
  x <- testset.gl
  set.seed(99)
  out <- gl.sample(x, nsample = 4, replace = TRUE, verbose = 0)
  src_pop <- as.character(pop(x))[match(src_names(out), indNames(x))]
  expect_identical(src_pop, as.character(pop(out)))
})

test_that("gl.sample default call covers every population", {
  x <- testset.gl
  # default nsample = min(table(pop(x))) = 1, because of the singleton pops
  expect_equal(min(table(pop(x))), 1)
  set.seed(7)
  out <- gl.sample(x, verbose = 0)
  singles <- names(which(table(pop(x)) == 1))
  # [approved F1]
  expect_true(all(singles %in% as.character(pop(out))))
  expect_equal(nInd(out), nPop(x))
})

# ---------------------------------------------------------------------------
# F2: @other survives the sampling and tracks the drawn individuals
# ---------------------------------------------------------------------------

test_that("gl.sample keeps @other and subsets it to the drawn individuals", {
  x <- platypus.gl
  expect_true(all(c("loc.metrics", "ind.metrics", "latlon",
                    "loc.metrics.flags", "history") %in% names(x@other)))

  set.seed(1)
  out <- gl.sample(x, nsample = 4, replace = TRUE, verbose = 0)

  # [approved F2] do.call(rbind, ...) dropped the whole @other list; a single
  # positive-index subset carries it
  expect_true(all(c("loc.metrics", "ind.metrics", "latlon",
                    "loc.metrics.flags", "history") %in% names(out@other)))
  expect_equal(nrow(out@other$ind.metrics), nInd(out))
  expect_equal(nrow(out@other$loc.metrics), nLoc(out))
  expect_equal(nrow(out@other$latlon), nInd(out))
  # every column of the input's individual metadata survives
  expect_equal(ncol(out@other$ind.metrics), ncol(x@other$ind.metrics))

  # rows track the drawn individuals exactly
  set.seed(1)
  samps <- ref_draw(x, 4, replace = TRUE)
  for (cn in setdiff(names(x@other$ind.metrics), "id")) {
    expect_equal(unname(out@other$ind.metrics[[cn]]),
                 unname(x@other$ind.metrics[samps, cn]), info = cn)
  }
  expect_equal(unname(as.matrix(out@other$latlon)),
               unname(as.matrix(x@other$latlon[samps, , drop = FALSE])))
  # loci are untouched, so loc.metrics is carried over verbatim
  expect_identical(out@other$loc.metrics, x@other$loc.metrics)
})

test_that("gl.sample tracks @other through duplicated draws", {
  x <- platypus.gl
  set.seed(4)
  out <- gl.sample(x, nsample = 20, replace = TRUE, onepop = TRUE, verbose = 0)
  set.seed(4)
  samps <- ref_draw(x, 20, replace = TRUE, onepop = TRUE)
  expect_gt(sum(duplicated(samps)), 0)

  # [approved F2] duplicated individuals get duplicated metadata rows
  expect_equal(nrow(out@other$ind.metrics), nInd(out))
  for (cn in setdiff(names(x@other$ind.metrics), "id")) {
    expect_equal(unname(out@other$ind.metrics[[cn]]),
                 unname(x@other$ind.metrics[samps, cn]), info = cn)
  }
  # duplicate source individuals do not produce duplicate names
  expect_false(any(duplicated(indNames(out))))
  expect_identical(src_names(out), indNames(x)[samps])
})

test_that("gl.sample output works with a standard downstream filter", {
  x <- platypus.gl
  set.seed(1)
  out <- gl.sample(x, nsample = 4, replace = TRUE, verbose = 0)
  # [approved F2] this raised "incorrect number of dimensions" while @other was
  # being destroyed
  f <- gl.filter.callrate(out, threshold = 0.5, plot.display = FALSE,
                          verbose = 0)
  expect_s4_class(f, "dartR")
  expect_equal(nInd(f), nInd(out))
  expect_lt(nLoc(f), nLoc(out))
})

test_that("gl.sample output needs no gl.compliance.check repair", {
  x <- platypus.gl
  set.seed(1)
  out <- gl.sample(x, nsample = 4, replace = TRUE, verbose = 0)
  fixed <- gl.compliance.check(out, verbose = 0)

  # [approved F2] previously gl.compliance.check was the documented workaround
  # and still could not recover ind.metrics beyond an id column, or latlon at
  # all; now there is nothing to recover
  expect_equal(ncol(out@other$ind.metrics), ncol(x@other$ind.metrics))
  expect_false(is.null(out@other$latlon))
  expect_equal(ncol(fixed@other$ind.metrics), ncol(x@other$ind.metrics))
  expect_false(is.null(fixed@other$latlon))
})

# ---------------------------------------------------------------------------
# F3: the two options(dartR_fbm) paths return the same contract
# ---------------------------------------------------------------------------

test_that("gl.sample returns the same metadata contract on both paths", {
  skip_if_not(requireNamespace("bigstatsr", quietly = TRUE),
              "bigstatsr not available for the FBM path")
  old <- getOption("dartR_fbm")
  on.exit(options(dartR_fbm = old), add = TRUE)

  x <- platypus.gl
  options(dartR_fbm = FALSE)
  set.seed(1)
  a <- gl.sample(x, nsample = 4, replace = TRUE, verbose = 0)

  options(dartR_fbm = TRUE)
  xf <- gl.gen2fbm(x, verbose = 0)
  set.seed(1)
  b <- gl.sample(xf, nsample = 4, replace = TRUE, verbose = 0)

  # [approved F3] identical call and seed under either global option
  expect_equal(sort(names(a@other)), sort(names(b@other)))
  expect_equal(nInd(a), nInd(b))
  expect_identical(indNames(a), indNames(b))
  expect_identical(as.character(pop(a)), as.character(pop(b)))
  expect_equal(dim(a@other$ind.metrics), dim(b@other$ind.metrics))
  expect_equal(unlist(a@other$loc.metrics.flags),
               unlist(b@other$loc.metrics.flags))
  expect_equal(length(a@other$history), length(b@other$history))
  expect_equal(unname(as.matrix(a)), unname(as.matrix(b)),
               ignore_attr = TRUE)
})

test_that("gl.sample keeps @other in sync on the FBM path", {
  skip_if_not(requireNamespace("bigstatsr", quietly = TRUE),
              "bigstatsr not available for the FBM path")
  old <- getOption("dartR_fbm")
  on.exit(options(dartR_fbm = old), add = TRUE)

  options(dartR_fbm = TRUE)
  xf <- gl.gen2fbm(platypus.gl, verbose = 0)
  set.seed(1)
  out <- gl.sample(xf, nsample = 4, replace = TRUE, verbose = 0)

  expect_true("ind.metrics" %in% names(out@other))
  expect_equal(nrow(out@other$ind.metrics), nInd(out))
})

# ---------------------------------------------------------------------------
# F4: locus metrics flags are honest after the individuals change
# ---------------------------------------------------------------------------

test_that("gl.sample flags locus metrics as stale", {
  skip_if_not(requireNamespace("bigstatsr", quietly = TRUE),
              "bigstatsr not available for the FBM path")
  old <- getOption("dartR_fbm")
  on.exit(options(dartR_fbm = old), add = TRUE)

  options(dartR_fbm = TRUE)
  xf <- gl.gen2fbm(testset.gl, verbose = 0)
  expect_true(xf@other$loc.metrics.flags$CallRate)

  set.seed(2)
  out <- gl.sample(xf, nsample = 6, replace = TRUE, onepop = TRUE, verbose = 0)

  # [approved F4] DAT4: the individuals changed, so every locus metric computed
  # across individuals is stale and every flag says so
  expect_false(out@other$loc.metrics.flags$CallRate)
  expect_true(all(unlist(out@other$loc.metrics.flags) == FALSE))

  # the stored values are deliberately left alone; the flags carry the warning
  recomputed <- 1 - colMeans(is.na(as.matrix(out)))
  stored <- out@other$loc.metrics$CallRate
  expect_gt(sum(abs(stored - recomputed) > 1e-9), 0)

  # so a downstream filter recalculates rather than trusting the stored values
  f <- gl.filter.callrate(out, threshold = 0.5, plot.display = FALSE,
                          verbose = 0)
  expect_true(f@other$loc.metrics.flags$CallRate)
})

test_that("gl.sample flags locus metrics as stale on the SNPbin path too", {
  x <- platypus.gl
  set.seed(1)
  out <- gl.sample(x, nsample = 4, replace = TRUE, verbose = 0)
  # [approved F4]
  expect_false(is.null(out@other$loc.metrics.flags))
  expect_true(all(unlist(out@other$loc.metrics.flags) == FALSE))
})

# ---------------------------------------------------------------------------
# F5: the call is appended to history, once, on both paths
# ---------------------------------------------------------------------------

test_that("gl.sample appends exactly one history entry", {
  x <- platypus.gl
  set.seed(1)
  out <- gl.sample(x, nsample = 4, replace = TRUE, verbose = 0)
  # [approved F5] FS8: the SNPbin path destroyed history outright and the FBM
  # path never extended it
  expect_length(out@other$history, length(x@other$history) + 1)
  expect_true(grepl(
    "^gl.sample",
    paste(deparse(out@other$history[[length(out@other$history)]]),
          collapse = " ")))

  skip_if_not(requireNamespace("bigstatsr", quietly = TRUE),
              "bigstatsr not available for the FBM path")
  old <- getOption("dartR_fbm")
  on.exit(options(dartR_fbm = old), add = TRUE)
  options(dartR_fbm = TRUE)
  xf <- gl.gen2fbm(x, verbose = 0)
  set.seed(1)
  outf <- gl.sample(xf, nsample = 4, replace = TRUE, verbose = 0)
  # [approved F5]
  expect_length(outf@other$history, length(xf@other$history) + 1)
})

# ---------------------------------------------------------------------------
# Individual naming
# ---------------------------------------------------------------------------

test_that("gl.sample makes individual names unique with a zero-padded prefix", {
  x <- platypus.gl
  set.seed(6)
  out <- gl.sample(x, nsample = 20, replace = TRUE, onepop = TRUE, verbose = 0)

  expect_false(any(duplicated(indNames(out))))
  expect_true(all(grepl("^[0-9]+_", indNames(out))))
  # prefix width matches the digit count of nInd
  expect_true(all(nchar(sub("_.*$", "", indNames(out))) ==
                    nchar(as.character(nInd(out)))))
})

test_that("gl.sample keeps ind.metrics$id in step with indNames", {
  skip_if_not(requireNamespace("bigstatsr", quietly = TRUE),
              "bigstatsr not available for the FBM path")
  old <- getOption("dartR_fbm")
  on.exit(options(dartR_fbm = old), add = TRUE)

  options(dartR_fbm = TRUE)
  xf <- gl.gen2fbm(testset.gl, verbose = 0)
  set.seed(2)
  out <- gl.sample(xf, nsample = 6, replace = TRUE, onepop = TRUE, verbose = 0)

  # [approved F6] indNames gained the NN_ prefix while ind.metrics$id did not
  expect_true(all(indNames(out) == as.character(out@other$ind.metrics$id)))
})

test_that("gl.sample keeps ind.metrics$id in step on the SNPbin path", {
  x <- platypus.gl
  set.seed(4)
  out <- gl.sample(x, nsample = 20, replace = TRUE, onepop = TRUE, verbose = 0)
  # [approved F6]
  expect_true(all(indNames(out) == as.character(out@other$ind.metrics$id)))
})

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

test_that("gl.sample handles SilicoDArT data and preserves ploidy 1", {
  x <- testset.gs
  set.seed(3)
  out <- gl.sample(x, nsample = 3, replace = TRUE, verbose = 0)

  expect_s4_class(out, "dartR")
  expect_true(all(ploidy(out) == 1))
  expect_equal(nLoc(out), nLoc(x))
  expect_true(all(as.matrix(out) %in% c(0, 1, NA)))
})

# ---------------------------------------------------------------------------
# Edge cases and argument validation
# ---------------------------------------------------------------------------

test_that("gl.sample nsample = 1 returns one individual per population", {
  x <- platypus.gl
  set.seed(8)
  out <- gl.sample(x, nsample = 1, replace = TRUE, verbose = 0)
  expect_equal(nInd(out), nPop(x))
})

test_that("gl.sample nsample equal to population size is allowed without replacement", {
  x <- platypus.gl
  smallest <- min(table(pop(x)))
  set.seed(9)
  out <- gl.sample(x, nsample = smallest, replace = FALSE, verbose = 0)
  expect_equal(nInd(out), smallest * nPop(x))
})

test_that("gl.sample rejects nsample = 0 by name", {
  x <- platypus.gl
  # [approved F7] previously "Subsetting resulted in zero individuals", which
  # names neither the argument nor the reason
  expect_error(gl.sample(x, nsample = 0, verbose = 0),
               "nsample must be 1 or more")
})

test_that("gl.sample rejects a proportion-style nsample", {
  x <- platypus.gl
  # [approved F7] a caller assuming nsample is a fraction is told it is a count
  expect_error(gl.sample(x, nsample = 0.5, verbose = 0),
               "nsample must be a whole number")
})

test_that("gl.sample rejects a fractional nsample instead of truncating it", {
  x <- platypus.gl
  # [approved F7] nsample = 2.7 silently returned 2 per population
  expect_error(gl.sample(x, nsample = 2.7, replace = TRUE, verbose = 0),
               "nsample must be a whole number")
})

test_that("gl.sample names the smallest population when nsample exceeds it", {
  x <- platypus.gl
  # [approved F7] previously the raw sample() error, naming no population
  expect_error(gl.sample(x, nsample = 25, replace = FALSE, verbose = 0),
               "smallest population")
  expect_error(gl.sample(x, nsample = 25, replace = FALSE, verbose = 0),
               names(which.min(table(pop(x)))))
})

test_that("gl.sample rejects a negative or NA nsample", {
  x <- platypus.gl
  expect_error(gl.sample(x, nsample = -1, verbose = 0))
  expect_error(gl.sample(x, nsample = NA, verbose = 0))
})

test_that("gl.sample rejects a non-logical replace or onepop", {
  x <- platypus.gl
  # [approved F7]
  expect_error(gl.sample(x, nsample = 2, replace = "yes", verbose = 0),
               "replace must be TRUE or FALSE")
  expect_error(gl.sample(x, nsample = 2, onepop = NA, verbose = 0),
               "onepop must be TRUE or FALSE")
})

test_that("gl.sample tolerates unused population factor levels", {
  x <- platypus.gl
  y <- x[pop(x) != "TENTERFIELD", ]
  set.seed(12)
  out <- gl.sample(y, nsample = 3, replace = TRUE, verbose = 0)
  expect_equal(nInd(out), 3 * length(unique(as.character(pop(y)))))
})

test_that("gl.sample tolerates an all-NA individual", {
  x <- platypus.gl[1:10, ]
  x@gen[[1]] <- new("SNPbin", rep(NA_integer_, nLoc(x)))
  set.seed(13)
  out <- gl.sample(x, nsample = 3, replace = TRUE, onepop = TRUE, verbose = 0)
  expect_equal(nInd(out), 3)
})

# ---------------------------------------------------------------------------
# Verbosity and purity
# ---------------------------------------------------------------------------

test_that("gl.sample is silent at verbose = 0", {
  x <- platypus.gl
  out <- capture.output(z <- gl.sample(x, nsample = 2, verbose = 0))
  expect_length(out, 0)
})

test_that("gl.sample reports completion at every verbosity above 0", {
  x <- platypus.gl
  for (v in c(1, 2, 3, 5)) {
    out <- capture.output(z <- gl.sample(x, nsample = 2, verbose = v))
    expect_true(any(grepl("Starting", out)), info = paste("verbose =", v))
    # [approved F8] FS9: no completion message printed at any verbosity
    expect_equal(sum(grepl("Completed", out)), 1L,
                 info = paste("verbose =", v))
  }
  # [approved F8] a progress line reports the draw from verbose 2
  out <- capture.output(z <- gl.sample(x, nsample = 2, verbose = 2))
  expect_true(any(grepl("Drew 2 individuals", out)))
  out <- capture.output(z <- gl.sample(x, nsample = 2, verbose = 1))
  expect_false(any(grepl("Drew", out)))
})

test_that("gl.sample does not mutate its input object", {
  x <- platypus.gl
  n_before <- nInd(x)
  pops_before <- as.character(pop(x))
  other_before <- names(x@other)

  set.seed(14)
  invisible(gl.sample(x, nsample = 4, replace = TRUE, onepop = TRUE, verbose = 0))

  expect_equal(nInd(x), n_before)
  expect_equal(as.character(pop(x)), pops_before)
  expect_equal(names(x@other), other_before)
})

test_that("gl.sample resolves the default nsample per onepop, as documented", {
  x <- platypus.gl
  set.seed(15); a <- gl.sample(x, verbose = 0)
  set.seed(15); b <- gl.sample(x, onepop = TRUE, verbose = 0)

  # onepop = FALSE: default is the smallest population, per population
  expect_equal(nInd(a), min(table(pop(x))) * nPop(x))
  # [approved F9] onepop = TRUE: default is nInd(x) as a whole-object total.
  # Unchanged in value, but now resolved deliberately rather than by the
  # evaluation order of a lazy promise against an overwritten pop(x).
  expect_equal(nInd(b), nInd(x))
})

test_that("gl.sample returns rows in draw order when nsample x nPop > nInd", {
  x <- platypus.gl
  nsample <- 30                     # 3 x 30 = 90 > 81
  set.seed(21)
  out <- gl.sample(x, nsample = nsample, replace = TRUE, verbose = 0)
  set.seed(21)
  samps <- ref_draw(x, nsample, replace = TRUE)

  expect_equal(nInd(out), length(samps))
  # [approved F10] the split/rbind assembly labelled chunks cyclically, so the
  # returned row order depended on whether nsample x nPop exceeded nInd
  expect_identical(as.character(pop(out)), as.character(pop(x))[samps])
  expect_identical(src_names(out), indNames(x)[samps])
  expect_identical(unname(as.matrix(out)),
                   unname(as.matrix(x)[samps, , drop = FALSE]))
  # metadata tracks the oversized draw too
  expect_equal(nrow(out@other$ind.metrics), nInd(out))
  expect_equal(nrow(out@other$loc.metrics), nLoc(out))
  expect_equal(unname(as.matrix(out@other$latlon)),
               unname(as.matrix(x@other$latlon[samps, , drop = FALSE])))

  # single-chunk case keeps draw order as it always did
  set.seed(21); out1 <- gl.sample(x, nsample = 5, replace = TRUE, verbose = 0)
  set.seed(21); s1 <- ref_draw(x, 5, replace = TRUE)
  expect_identical(as.character(pop(out1)), as.character(pop(x))[s1])
})

# ---------------------------------------------------------------------------
# Negative control: the SNPbin duplicate-index defect does NOT apply here
# ---------------------------------------------------------------------------

test_that("resampling individuals with replacement preserves NA integrity", {
  # The known adegenet SNPbin[] duplicate-index defect corrupts NA when LOCI
  # (columns) are subset with repeated indices. gl.sample duplicates
  # INDIVIDUALS (rows), which selects whole SNPbin objects from @gen and never
  # calls SNPbin's own [ method. This test guards that distinction.
  x <- testset.gl
  idx <- c(5, 5, 5, 7, 7, 9)
  sub <- x[idx, ]
  expect_equal(sum(is.na(as.matrix(sub))),
               sum(is.na(as.matrix(x)[idx, , drop = FALSE])))
  expect_identical(unname(as.matrix(sub)),
                   unname(as.matrix(x)[idx, , drop = FALSE]))

  set.seed(1234)
  out <- gl.sample(x, nsample = 5, replace = TRUE, verbose = 0)
  # [approved F1] the reference draw uses the fixed idiom, so it selects the
  # rows gl.sample now selects
  set.seed(1234)
  samps <- ref_draw(x, 5, replace = TRUE)
  expect_equal(sum(is.na(as.matrix(out))),
               sum(is.na(as.matrix(x)[samps, , drop = FALSE])))
  expect_identical(unname(as.matrix(out)),
                   unname(as.matrix(x)[samps, , drop = FALSE]))
})
