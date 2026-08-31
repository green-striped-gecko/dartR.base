# Characterization tests for gl.make.recode.pop
# Baseline snapshotted before review (review-gl.make.recode.pop).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.make.recode.pop writes an exact two-column proforma", {
  tmpd <- tempdir()
  f <- file.path(tmpd, "chartest_rp.csv")
  out <- capture.output(
    res <- gl.make.recode.pop(testset.gl, out.recode.file = "chartest_rp.csv",
                              outpath = tmpd, verbose = 0)
  )
  expect_length(out, 0)
  tab <- read.csv(f, header = FALSE, stringsAsFactors = FALSE)
  expect_equal(nrow(tab), nPop(testset.gl))
  expect_identical(tab[, 1], tab[, 2])
  expect_setequal(tab[, 1], levels(pop(testset.gl)))
  unlink(f)
})

test_that("gl.make.recode.pop proforma round-trips through gl.recode.pop", {
  tmpd <- tempdir()
  f <- file.path(tmpd, "chartest_rp_rt.csv")
  invisible(gl.make.recode.pop(testset.gl, out.recode.file = "chartest_rp_rt.csv",
                               outpath = tmpd, verbose = 0))
  res <- gl.recode.pop(testset.gl, pop.recode = f, verbose = 0)
  expect_identical(as.character(pop(res)), as.character(pop(testset.gl)))
  unlink(f)
})

test_that("gl.make.recode.pop error path and SilicoDArT", {
  xnp <- testset.gl
  xnp@pop <- NULL
  expect_error(
    gl.make.recode.pop(xnp, outpath = tempdir(), verbose = 0),
    "Population names not detected"
  )
  tmpd <- tempdir()
  f <- file.path(tmpd, "chartest_rp_gs.csv")
  xcopy <- testset.gs
  invisible(gl.make.recode.pop(xcopy, out.recode.file = "chartest_rp_gs.csv",
                               outpath = tmpd, verbose = 0))
  expect_identical(xcopy, testset.gs)
  expect_equal(nrow(read.csv(f, header = FALSE)), nPop(testset.gs))
  unlink(f)
})

test_that("gl.make.recode.pop return value", {
  # [approved diff F1] returns NULL invisibly; docs state the file-writing
  # contract.
  tmpd <- tempdir()
  v <- withVisible(gl.make.recode.pop(testset.gl,
                                      out.recode.file = "chartest_rp_v.csv",
                                      outpath = tmpd, verbose = 0))
  expect_null(v$value)
  expect_false(v$visible)
  unlink(file.path(tmpd, "chartest_rp_v.csv"))
})
