# Characterization tests for gl.merge.pop
# Baseline snapshotted before review (review-gl.merge.pop).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("two populations merge into one", {
  o <- capture.output(g <- gl.merge.pop(testset.gl,
        old = c("EmsubRopeMata", "EmvicVictJasp"), new = "Outgroup",
        verbose = 0))
  expect_length(o, 0)
  expect_true("Outgroup" %in% popNames(g))
  expect_equal(nPop(g), nPop(testset.gl) - 1)
  expect_equal(sum(pop(g) == "Outgroup"),
               sum(pop(testset.gl) %in% c("EmsubRopeMata", "EmvicVictJasp")))
  expect_equal(length(g@other$history),
               length(testset.gl@other$history) + 1)
})

test_that("single old population is a rename", {
  invisible(capture.output(g <- gl.merge.pop(testset.gl,
        old = "EmsubRopeMata", new = "Renamed", verbose = 0)))
  expect_true("Renamed" %in% popNames(g))
  expect_false("EmsubRopeMata" %in% popNames(g))
  expect_equal(nPop(g), nPop(testset.gl))
})

test_that("nonexistent old population", {
  # [approved diff F2] baseline: silently ignored - no error, no
  # message, object unchanged. Now fatal, naming the missing pops.
  expect_error(capture.output(gl.merge.pop(testset.gl, old = "NOSUCHPOP",
                                           new = "Merged", verbose = 0)),
               "not present")
})

test_that("empty and NULL old are fatal at every verbosity", {
  expect_error(capture.output(gl.merge.pop(testset.gl, old = character(0),
                                           new = "M", verbose = 2)))
  # [approved diff F1] baseline: at verbose 0 the empty-old validation
  # sat inside the verbose >= 1 gate and the call silently no-opped.
  expect_error(capture.output(gl.merge.pop(testset.gl,
        old = character(0), new = "M", verbose = 0)))
  expect_error(capture.output(gl.merge.pop(testset.gl, old = NULL,
                                           new = "M", verbose = 0)))
  expect_error(capture.output(gl.merge.pop(testset.gl,
        old = "EmsubRopeMata", new = NULL, verbose = 0)))
})
