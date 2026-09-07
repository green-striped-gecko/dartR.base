# Characterization tests for gl.define.pop
# Baseline snapshotted before review (review-gl.define.pop).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("individuals assigned to the new population", {
  o <- capture.output(g <- gl.define.pop(testset.gl,
        ind.list = c("AA019073", "AA004859"), new = "newguys", verbose = 0))
  expect_length(o, 0)
  expect_equal(sum(pop(g) == "newguys"), 2)
  expect_true(all(indNames(g)[pop(g) == "newguys"] %in%
                    c("AA019073", "AA004859")))
  expect_equal(nInd(g), nInd(testset.gl))
  expect_false("newguys" %in% popNames(testset.gl))   # input untouched
  # single history append, stored as a call
  expect_equal(length(g@other$history),
               length(testset.gl@other$history) + 1)
  expect_true(is.call(g@other$history[[length(g@other$history)]]))
})

test_that("not-present individuals are ignored with a gated warning", {
  # [approved diff F1] baseline: the warning printed at verbose 0.
  o <- capture.output(g <- gl.define.pop(testset.gl,
        ind.list = c("AA019073", "NOSUCH"), new = "newguys", verbose = 0))
  expect_length(o, 0)
  expect_equal(sum(pop(g) == "newguys"), 1)
})

test_that("all-absent ind.list is fatal", {
  expect_error(suppressWarnings(capture.output(
    gl.define.pop(testset.gl, ind.list = "NOSUCH", new = "x", verbose = 0))))
})
