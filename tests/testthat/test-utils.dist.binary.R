# Characterization tests for utils.dist.binary
# Baseline snapshotted before review (distance-wave review, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.

test_that("jaccard and sorensen match manual N-count recomputations", {
  sub <- gl.keep.ind(testset.gs, ind.list = indNames(testset.gs)[1:5],
                     verbose = 0)
  sm <- as.matrix(sub)
  o <- capture.output({
    dj <- as.matrix(utils.dist.binary(sub, method = "jaccard", verbose = 0))
    ds <- as.matrix(utils.dist.binary(sub, method = "sorensen", verbose = 0))
  })
  v1 <- sm[1, ]; v2 <- sm[2, ]; ok <- !is.na(v1) & !is.na(v2)
  N11 <- sum(v1[ok] == 1 & v2[ok] == 1)
  N10 <- sum(v1[ok] == 1 & v2[ok] == 0)
  N01 <- sum(v1[ok] == 0 & v2[ok] == 1)
  expect_equal(dj[1, 2], (N01 + N10) / (N11 + N01 + N10))
  expect_equal(ds[1, 2], (N01 + N10) / (2 * N11 + N01 + N10))
})

test_that("bray-curtis is accepted as documented", {
  sub <- gl.keep.ind(testset.gs, ind.list = indNames(testset.gs)[1:5],
                     verbose = 0)
  # [approved diff D5] baseline: bray-curtis was documented and implemented
  # but missing from the validation list, so it silently fell back to
  # simple matching.
  o <- capture.output(db <- utils.dist.binary(sub, method = "bray-curtis",
                                              verbose = 2))
  expect_false(any(grepl("not in the list", o)))  # [approved diff D5]
})

test_that("the scale warning respects verbose 0", {
  sub <- gl.keep.ind(testset.gs, ind.list = indNames(testset.gs)[1:5],
                     verbose = 0)
  # [approved diff D6] baseline: leaked at verbose 0.
  o <- capture.output(d <- utils.dist.binary(sub, method = "jaccard",
                                             scale = TRUE, verbose = 0))
  expect_equal(length(o), 0)  # [approved diff D6]
})
