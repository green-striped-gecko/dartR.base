test_that("gl.keep.loc retains only the nominated loci (SNP)", {
  x <- testset.gl
  to_keep <- locNames(x)[c(3, 7)]
  x2 <- gl.keep.loc(x, loc.list = to_keep, verbose = 0)
  expect_setequal(locNames(x2), to_keep)
  expect_equal(nInd(x2), nInd(x))
  expect_true(all(ploidy(x2) == 2))
  expect_length(x2@other$history, length(x@other$history) + 1)
})

test_that("gl.keep.loc works on SilicoDArT data", {
  x <- testset.gs
  to_keep <- locNames(x)[1:2]
  x2 <- gl.keep.loc(x, loc.list = to_keep, verbose = 0)
  expect_setequal(locNames(x2), to_keep)
})

test_that("gl.keep.loc ignores a nonexistent locus name with a warning, not fatal", {
  x <- testset.gl
  real_loc <- locNames(x)[1]
  expect_error(
    x2 <- gl.keep.loc(x, loc.list = c(real_loc, "NOT_A_REAL_LOCUS"), verbose = 0),
    NA
  )
  expect_true(real_loc %in% locNames(x2))
})

test_that("gl.keep.loc: first/last range within bounds works", {
  x <- testset.gl
  x2 <- gl.keep.loc(x, first = 2, last = 5, verbose = 0)
  expect_equal(nLoc(x2), 4)
  expect_setequal(locNames(x2), locNames(x)[2:5])
})

test_that("neither loc.list nor first specified errors clearly (F2 fix)", {
  x <- testset.gl
  expect_error(
    gl.keep.loc(x, verbose = 0),
    "Need to specify either a range of loci"
  )
})

test_that("first given, last omitted defaults to the last locus (F1 fix)", {
  x <- testset.gl
  n <- nLoc(x)
  x2 <- gl.keep.loc(x, first = 5, verbose = 0)
  expect_equal(nLoc(x2), n - 4)
  expect_setequal(locNames(x2), locNames(x)[5:n])
})

test_that("last > nLoc(x) warns and clamps; first > nLoc(x) is fatal (F3 fix); clamp warnings silent at verbose=0", {
  x <- testset.gl
  n <- nLoc(x)
  out <- capture.output(x2 <- gl.keep.loc(x, first = n - 2, last = n + 500, verbose = 0))
  expect_length(out, 0)
  expect_equal(nLoc(x2), 3)
  expect_setequal(locNames(x2), locNames(x)[(n - 2):n])

  expect_error(
    gl.keep.loc(x, first = n + 50, last = n + 200, verbose = 0),
    "cannot be greater than the number of loci"
  )
})

test_that("zero-length loc.list still returns the object unchanged (external callers depend on this)", {
  # dartR.popgen (gl.ld.haplotype) and dartR.sim (gl.sim.WF.run) pass
  # character(0) loc.lists and rely on the warn-and-return-unchanged path;
  # the F2 fatal error must fire only when loc.list is NULL, not empty.
  x <- testset.gl
  expect_error(x2 <- gl.keep.loc(x, loc.list = character(0), verbose = 0), NA)
  expect_equal(nLoc(x2), nLoc(x))
})
