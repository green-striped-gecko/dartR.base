# Characterization tests for utils.jackknife
# Baseline snapshotted before review (distance-wave review, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.

make_g <- function() {
  set.seed(21)
  m <- matrix(sample(0:2, 10 * 30, replace = TRUE), nrow = 10)
  g <- new("genlight", gen = m, ploidy = 2)
  indNames(g) <- paste0("i", 1:10); locNames(g) <- paste0("L", 1:30)
  pop(g) <- factor(rep(c("A", "B"), each = 5))
  suppressWarnings(gl.compliance.check(g, verbose = 0))
}

test_that("pop jackknife returns one result per unit, silently at verbose 0", {
  g <- make_g()
  # [approved diff D11] baseline: the gl.set.verbosity save/restore inside
  # each replicate leaked 6 lines at verbose 0.
  o <- capture.output(jk <- utils.jackknife(g, FUN = "gl.alf", unit = "pop",
        recalc = FALSE, mono.rm = FALSE, n.cores = 1, verbose = 0))
  expect_equal(length(jk), 2)
  expect_equal(length(o), 0)  # [approved diff D11]
  expect_true(is.data.frame(jk[[1]]))
})

test_that("a unit vector of length > 1 fails informatively", {
  g <- make_g()
  # [approved diff D10] baseline: length(unit == 1) - the length of the
  # comparison, not of unit - crashed with "the condition has length > 1"
  # instead of reaching the informative stop.
  expect_error(
    utils.jackknife(g, FUN = "gl.alf", unit = c("loc", "ind"), verbose = 0),
    "length 1")  # [approved diff D10]
})

test_that("an invalid unit fails informatively", {
  g <- make_g()
  expect_error(
    utils.jackknife(g, FUN = "gl.alf", unit = "chromosome", verbose = 0),
    "loc")
})
