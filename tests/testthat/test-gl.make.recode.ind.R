# Characterization tests for gl.make.recode.ind
# Baseline snapshotted before review (review-gl.make.recode.ind).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.make.recode.ind writes an exact two-column proforma", {
  tmpd <- tempdir()
  f <- file.path(tmpd, "chartest_ri.csv")
  out <- capture.output(
    res <- gl.make.recode.ind(testset.gl, out.recode.file = "chartest_ri.csv",
                              outpath = tmpd, verbose = 0)
  )
  expect_length(out, 0)
  expect_true(file.exists(f))
  tab <- read.csv(f, header = FALSE, stringsAsFactors = FALSE)
  expect_equal(nrow(tab), nInd(testset.gl))
  expect_identical(tab[, 1], tab[, 2])
  expect_identical(tab[, 1], indNames(testset.gl))
  unlink(f)
})

test_that("gl.make.recode.ind proforma round-trips through gl.recode.ind", {
  tmpd <- tempdir()
  f <- file.path(tmpd, "chartest_ri_rt.csv")
  invisible(gl.make.recode.ind(testset.gl, out.recode.file = "chartest_ri_rt.csv",
                               outpath = tmpd, verbose = 0))
  res <- gl.recode.ind(testset.gl, ind.recode = f, verbose = 0)
  expect_identical(indNames(res), indNames(testset.gl))
  unlink(f)
})

test_that("gl.make.recode.ind works on SilicoDArT and leaves input untouched", {
  tmpd <- tempdir()
  f <- file.path(tmpd, "chartest_ri_gs.csv")
  xcopy <- testset.gs
  invisible(gl.make.recode.ind(xcopy, out.recode.file = "chartest_ri_gs.csv",
                               outpath = tmpd, verbose = 0))
  expect_identical(xcopy, testset.gs)
  tab <- read.csv(f, header = FALSE)
  expect_equal(nrow(tab), nInd(testset.gs))
  unlink(f)
})

test_that("gl.make.recode.ind return value", {
  # [approved diff F1] returns NULL invisibly (no 'NULL' auto-print);
  # docs now state the file-writing contract.
  tmpd <- tempdir()
  v <- withVisible(gl.make.recode.ind(testset.gl,
                                      out.recode.file = "chartest_ri_v.csv",
                                      outpath = tmpd, verbose = 0))
  expect_null(v$value)
  expect_false(v$visible)
  unlink(file.path(tmpd, "chartest_ri_v.csv"))
})
