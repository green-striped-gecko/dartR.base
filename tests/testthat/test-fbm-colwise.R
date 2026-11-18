# Column-wise statistics of FBM-backed objects are computed in blocks of
# loci (.gl_colwise) by gl.compliance.check, the utils.recalc.* helpers and
# gl.alf. The values must be identical to those of the gen-backed object,
# for any block size, and the whole matrix must never be decoded.

fbm_pair <- function() {
  x <- gl.compliance.check(testset.gl[1:30, 1:60], verbose = 0)
  m <- as.matrix(x)
  m[, 5] <- NA                               # an all-NA locus
  x@gen <- new("genlight", m, ploidy = ploidy(x))@gen
  list(x = x, f = suppressWarnings(gl.gen2fbm(x, verbose = 0)))
}

test_that("locus metrics and gl.alf are identical for gen and FBM objects", {
  skip_if_not_installed("bigstatsr")
  old <- options(bigstatsr.block.sizeGB = getOption("bigstatsr.block.sizeGB"))
  on.exit(options(old))
  xs <- fbm_pair()
  for (bs in c(1, 1e-6)) {                   # 1e-6: one locus per block
    options(bigstatsr.block.sizeGB = bs)
    expect_identical(gl.alf(xs$f, verbose = 0), gl.alf(xs$x, verbose = 0))
    expect_identical(gl.recalc.metrics(xs$f, verbose = 0)@other$loc.metrics,
                     gl.recalc.metrics(xs$x, verbose = 0)@other$loc.metrics)
    cf <- gl.compliance.check(xs$f, verbose = 0)
    cg <- gl.compliance.check(xs$x, verbose = 0)
    expect_identical(cf@other$loc.metrics, cg@other$loc.metrics)
    expect_identical(cf@other$loc.metrics.flags, cg@other$loc.metrics.flags)
    expect_identical(.gl_colwise(xs$f, function(t) colSums(is.na(t))),
                     colSums(is.na(as.matrix(xs$x))))
  }
})

test_that("the value check of gl.compliance.check reads FBM blocks", {
  skip_if_not_installed("bigstatsr")
  xs <- fbm_pair()
  expect_output(gl.compliance.check(xs$f, verbose = 1), "0, 1 or 2 confirmed")
  # a code in which raw 3 is NA (as in CODE_012) and raw 5 decodes to 5
  bad <- xs$f
  G <- bigstatsr::FBM.code256(nInd(bad), nLoc(bad),
                              code = c(0, 1, 2, NA, 4, 5, rep(NA, 250)),
                              backingfile = tempfile())
  m <- bad@fbm[]
  m[is.na(m)] <- 3
  m[1, 1] <- 5
  G[] <- m
  bad@fbm <- G
  expect_output(gl.compliance.check(bad, verbose = 1), "must be scored")
})

test_that("gl.impute and gl.compliance.check never decode the whole FBM", {
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  xs <- fbm_pair()
  local_mocked_bindings(.fbm_as_matrix = function(x) {
    stop("the whole FBM was decoded")
  })
  expect_no_error(gl.compliance.check(xs$f, verbose = 0))
  for (meth in c("frequency", "HW", "neighbour", "random")) {
    set.seed(1)
    expect_no_error(capture.output(
      gl.impute(xs$f, method = meth, verbose = 3)))
  }
})
