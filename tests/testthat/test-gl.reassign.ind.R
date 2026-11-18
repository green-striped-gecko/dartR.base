# Characterization tests for gl.reassign.ind
# Baseline snapshotted before review (review-gl.reassign.ind).

test_that("character, numeric and logical selections reassign correctly", {
  o <- capture.output(g <- gl.reassign.ind(testset.gl,
        ind.list = indNames(testset.gl)[1:5], new.pop = "NewGroup",
        verbose = 0))
  expect_length(o, 0)
  expect_equal(sum(pop(g) == "NewGroup"), 5)
  invisible(capture.output(g2 <- gl.reassign.ind(testset.gl,
        ind.list = 1:10, new.pop = "Grp", verbose = 0)))
  expect_equal(sum(pop(g2) == "Grp"), 10)
  sel <- rep(FALSE, nInd(testset.gl)); sel[c(3, 7)] <- TRUE
  invisible(capture.output(g3 <- gl.reassign.ind(testset.gl,
        ind.list = sel, new.pop = "Two", verbose = 0)))
  expect_equal(sum(pop(g3) == "Two"), 2)
  # others retain their assignment
  expect_equal(as.character(pop(g3))[!sel], as.character(pop(testset.gl))[!sel])
  expect_equal(length(g@other$history),
               length(testset.gl@other$history) + 1)
})

test_that("invalid selections are fatal", {
  expect_error(gl.reassign.ind(testset.gl, ind.list = "NOSUCH",
                               new.pop = "X", verbose = 0))
  expect_error(gl.reassign.ind(testset.gl, ind.list = c(0, 5),
                               new.pop = "X", verbose = 0))
  expect_error(gl.reassign.ind(testset.gl, ind.list = c(TRUE, FALSE),
                               new.pop = "X", verbose = 0))
  expect_error(gl.reassign.ind(testset.gl, ind.list = 1:3,
                               new.pop = "", verbose = 0))
})

test_that("reassignment to an existing population works", {
  target <- popNames(testset.gl)[2]
  before <- sum(pop(testset.gl) == target)
  invisible(capture.output(g <- gl.reassign.ind(testset.gl, ind.list = 1:3,
        new.pop = target, verbose = 0)))
  expect_gte(sum(pop(g) == target), before)
  expect_equal(nInd(g), nInd(testset.gl))
})
