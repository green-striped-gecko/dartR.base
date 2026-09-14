# Characterization tests for utils.n.var.invariant
# Baseline snapshotted before review (kernel-wave review, dev at ddaed27).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("variant/invariant site counts verify on platypus", {
  skip_if_not_installed("dartR.data")
  library(dartR.data)
  o <- capture.output(nv <- utils.n.var.invariant(platypus.gl, verbose = 0))
  expect_equal(length(o), 0)
  lm <- nv@other$loc.metrics
  tab <- table(lm$CloneID)
  expect_true(all(lm$n.variant == as.integer(tab[as.character(lm$CloneID)])))
  expect_true(all(lm$n.invariant == lm$lenTrimSeq - lm$n.variant, na.rm = TRUE))
  h <- nv@other$history
  expect_true(is.call(h[[length(h)]]))
})

test_that("the secondaries-history warning respects verbose 0", {
  skip_if_not_installed("dartR.data")
  library(dartR.data)
  pg <- platypus.gl
  pg@other$history <- list(quote(gl.filter.secondaries(x)))
  # [approved diff K8] baseline: the warning printed at verbose 0.
  o <- capture.output(nv <- utils.n.var.invariant(pg, verbose = 0))
  expect_equal(length(o), 0)  # [approved diff K8]
})
