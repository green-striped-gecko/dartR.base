# Characterization tests for utils.collapse.matrix
# Baseline snapshotted before review (infrastructure wave, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.

test_that("off-diagonal population means are exact; diagonal excludes self", {
  sub <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1:3],
                     verbose = 0)
  o <- capture.output({
    Dind <- gl.dist.ind(sub, method = "euclidean", verbose = 0)
    Dpop <- utils.collapse.matrix(D = as.matrix(Dind), x = sub, verbose = 0)
  })
  mat <- as.matrix(Dind)
  i1 <- indNames(sub)[pop(sub) == popNames(sub)[1]]
  i2 <- indNames(sub)[pop(sub) == popNames(sub)[2]]
  expect_equal(unname(Dpop[2, 1]), mean(mat[i2, i1]))
  # [approved diff I9] baseline: the within-population mean included the
  # zero self-distances (1.6195 including vs 1.7995 excluding on this
  # fixture).
  m11 <- mat[i1, i1]
  expect_equal(unname(Dpop[1, 1]), mean(m11[upper.tri(m11)]))  # [approved diff I9]
})

test_that("verbose 0 is silent and dist input returns dist", {
  sub <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1:3],
                     verbose = 0)
  o0 <- capture.output(Dind <- gl.dist.ind(sub, method = "euclidean",
                                           verbose = 0))
  o <- capture.output(Dpop <- utils.collapse.matrix(D = Dind, x = sub,
                                                    verbose = 0))
  expect_equal(length(o), 0)
  expect_s3_class(Dpop, "dist")
})
