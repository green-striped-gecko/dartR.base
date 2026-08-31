test_that("gl.drop.loc removes the nominated loci (SNP)", {
  x <- testset.gl
  to_drop <- locNames(x)[c(3, 7)]
  x2 <- gl.drop.loc(x, loc.list = to_drop, verbose = 0)
  expect_false(any(to_drop %in% locNames(x2)))
  expect_equal(nLoc(x2), nLoc(x) - 2)
  expect_equal(nInd(x2), nInd(x))
  expect_true(all(ploidy(x2) == 2))
  expect_length(x2@other$history, length(x@other$history) + 1)
  expect_equal(nrow(x2@other$loc.metrics), nLoc(x2))
})

test_that("gl.drop.loc works on SilicoDArT data", {
  x <- testset.gs
  to_drop <- locNames(x)[1:2]
  x2 <- gl.drop.loc(x, loc.list = to_drop, verbose = 0)
  expect_false(any(to_drop %in% locNames(x2)))
  expect_equal(nLoc(x2), nLoc(x) - 2)
})

test_that("gl.drop.loc: range within bounds drops that range", {
  x <- testset.gl
  x2 <- gl.drop.loc(x, first = 2, last = 5, verbose = 0)
  expect_equal(nLoc(x2), nLoc(x) - 4)
  expect_false(any(locNames(x)[2:5] %in% locNames(x2)))
})

test_that("gl.drop.loc errors clearly when neither loc.list nor first is given", {
  x <- testset.gl
  expect_error(
    gl.drop.loc(x, verbose = 0),
    "Need to specify either a range of loci"
  )
})

test_that("gl.drop.loc: nonexistent locus in loc.list is ignored, real one still dropped", {
  x <- testset.gl
  real_loc <- locNames(x)[1]
  expect_error(
    x2 <- gl.drop.loc(x, loc.list = c(real_loc, "NOT_A_REAL_LOCUS"), verbose = 0),
    NA
  )
  expect_false(real_loc %in% locNames(x2))
  expect_equal(nLoc(x2), nLoc(x) - 1)
})

test_that("the not-present warning names the loci the user listed, not dataset loci (F1 fix)", {
  x <- testset.gl
  out <- capture.output(
    x2 <- gl.drop.loc(x, loc.list = c("BOGUS_LOCUS_XYZ", locNames(x)[5]), verbose = 2)
  )
  warn_line <- grep("not present", out, value = TRUE)
  expect_length(warn_line, 1)
  expect_true(grepl("BOGUS_LOCUS_XYZ", warn_line))
  expect_false(grepl(locNames(x)[1], warn_line, fixed = TRUE))
  expect_equal(nLoc(x2), nLoc(x) - 1)
})

test_that("first given, last omitted defaults to the last locus (F2 fix)", {
  x <- testset.gl
  n <- nLoc(x)
  x2 <- gl.drop.loc(x, first = 250, verbose = 0)
  expect_equal(nLoc(x2), 249)
  expect_setequal(locNames(x2), locNames(x)[1:249])
})

test_that("last > nLoc warns and clamps; first > nLoc is fatal (F3 fix); clamp warnings silent at verbose=0 (F7)", {
  x <- testset.gl
  n <- nLoc(x)
  out <- capture.output(x2 <- gl.drop.loc(x, first = n - 2, last = n + 500, verbose = 0))
  expect_length(out, 0)
  expect_equal(nLoc(x2), n - 3)
  expect_error(
    gl.drop.loc(x, first = n + 50, last = n + 200, verbose = 0),
    "cannot be greater than the number of loci"
  )
})
