# Characterization tests for gl.filter.replicates
# Baseline snapshotted before review (review-gl.filter.replicates).
# Assertions marked [approved diff] were flipped in Phase C.

skip_if_not_installed("Rcpp")
skip_if_not_installed("RcppParallel")

make_dup2 <- function() {
  base <- platypus.gl[1:30, 1:500]
  bm <- as.matrix(base)
  dm <- rbind(bm, bm[1:2, ])
  rownames(dm) <- make.unique(c(rownames(bm), rownames(bm)[1:2]))
  dup <- new("genlight", gen = dm, ploidy = 2)
  indNames(dup) <- rownames(dm)
  locNames(dup) <- colnames(bm)
  pop(dup) <- factor(rep("A", nrow(dm)))
  suppressWarnings(gl.compliance.check(dup, verbose = 0))
}
dup <- make_dup2()
invisible(capture.output(rep.res <- gl.report.replicates(dup,
      loc_threshold = 100, perc_geno = 0.9, plot.out = FALSE,
      verbose = 0)))

test_that("filter removes one member per replicate pair", {
  invisible(capture.output(f <- gl.filter.replicates(dup,
        replicates.report = rep.res, loc_threshold = 100,
        perc_geno = 0.9, verbose = 0)))
  # [approved diff F1] baseline: tied-missingness pairs lost BOTH
  # members - 5 individuals dropped where 3 (one per pair) is correct.
  expect_equal(nInd(f), 29)
  expect_true("T27" %in% indNames(f))     # tie: earlier name kept
  expect_false("T27.1" %in% indNames(f))
  expect_true("T35" %in% indNames(f))
  expect_false("T35.1" %in% indNames(f))
  expect_false("T3" %in% indNames(f))     # correct: more missing than T5
  expect_true("T5" %in% indNames(f))
  # [approved diff F3] baseline: the last history entry was
  # gl.drop.ind's call.
  h <- f@other$history[[length(f@other$history)]]
  expect_equal(deparse(h[[1]]), "gl.filter.replicates")
})

test_that("string report from the no-pairs path", {
  # [approved diff F2] baseline: crashed with "$ operator is invalid
  # for atomic vectors". Now a clear validation error.
  expect_error(capture.output(gl.filter.replicates(dup,
        replicates.report = "No pair of individuals ...", verbose = 0)),
        "gl.report.replicates")
})

test_that("re-thresholding to an empty set", {
  # [approved diff F1] baseline: crashed via gl.drop.ind ("no
  # individuals to drop"). Now returns the object unchanged.
  o <- capture.output(f <- gl.filter.replicates(dup,
        replicates.report = rep.res, loc_threshold = 499,
        perc_geno = 0.9999999, verbose = 0))
  expect_length(o, 0)
  expect_equal(nInd(f), nInd(dup))
})
