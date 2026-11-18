# Characterization tests for glSum and glMean (dartR overrides of the
# adegenet functions), captured before the function-review changes.

x <- testset.gl

test_that("dense objects delegate to adegenet", {
  for (a in c(TRUE, FALSE)) {
    g <- x
    class(g) <- "genlight"
    expect_identical(glSum(x, alleleAsUnit = a),
                     adegenet::glSum(g, alleleAsUnit = a))
    expect_identical(glMean(x, alleleAsUnit = a),
                     adegenet::glMean(g, alleleAsUnit = a))
  }
  expect_true(is(x, "dartR"))
})

test_that("glMean equals an independent allele frequency", {
  m <- as.matrix(x)
  q <- colSums(m, na.rm = TRUE) / (2 * colSums(!is.na(m)))
  expect_equal(unname(glMean(x)), unname(q))
  expect_equal(unname(glSum(x)), unname(colSums(m, na.rm = TRUE)))
  expect_equal(sum(is.nan(glMean(x[1:5, ]))), 15L)
})

test_that("FBM-backed objects give the dense result", {
  skip_if_not_installed("bigstatsr")
  xf <- gl.gen2fbm(x, verbose = 0)
  for (idx in list(list(i = seq_len(nInd(x)), j = seq_len(nLoc(x))),
                   list(i = 1:20, j = 1:50),
                   list(i = nInd(x):1, j = seq_len(nLoc(x))),
                   list(i = 1:5, j = seq_len(nLoc(x))))) {
    d <- x[idx$i, idx$j]
    f <- xf[idx$i, idx$j]
    for (a in c(TRUE, FALSE)) {
      expect_equal(glSum(f, alleleAsUnit = a), glSum(d, alleleAsUnit = a))
      expect_identical(typeof(glSum(f, alleleAsUnit = a)),
                       typeof(glSum(d, alleleAsUnit = a)))
      expect_equal(glMean(f, alleleAsUnit = a), glMean(d, alleleAsUnit = a))
    }
  }
})

test_that("plain genlight input", {
  g <- new("genlight", as.matrix(x)[1:10, 1:20], ploidy = 2)
  expect_identical(glMean(g), adegenet::glMean(g))
  expect_identical(glSum(g), adegenet::glSum(g))
})
