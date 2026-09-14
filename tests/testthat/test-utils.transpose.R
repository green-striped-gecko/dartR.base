# Characterization tests for utils.transpose
# Baseline snapshotted before review (infrastructure wave, dev at ddaed27).

test_that("transpose swaps dimensions, names and metrics; round-trip exact", {
  sub <- gl.compliance.check(testset.gl[1:8, 1:20], verbose = 0)
  tr <- utils.transpose(sub)
  expect_equal(nInd(tr), nLoc(sub))
  expect_equal(nLoc(tr), nInd(sub))
  expect_equal(indNames(tr), locNames(sub))
  expect_equal(locNames(tr), indNames(sub))
  expect_equal(as.matrix(utils.transpose(tr)), as.matrix(sub))
})
