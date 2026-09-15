# Characterization tests for gl.write.csv
# Baseline snapshotted before review (review-gl.write.csv).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.write.csv produces the documented structure for SNP data", {
  td <- tempdir()
  f <- file.path(td, "wcsv-snp.csv")
  out <- capture.output(
    gl.write.csv(testset.gl, outfile = "wcsv-snp.csv", outpath = td,
                 verbose = 0)
  )
  raw <- readLines(f)
  expect_equal(length(raw), nLoc(testset.gl) + 2)
  hdr <- gsub('"', '', strsplit(raw[1], ",")[[1]])
  row1 <- gsub('"', '', strsplit(raw[2], ",")[[1]])
  nmeta <- ncol(testset.gl@other$loc.metrics)
  expect_equal(length(hdr), nmeta + nInd(testset.gl))
  expect_equal(hdr[1:nmeta], colnames(testset.gl@other$loc.metrics))
  expect_equal(hdr[nmeta + 1], indNames(testset.gl)[1])
  expect_true(all(row1[1:nmeta] == "*"))
  expect_equal(row1[nmeta + 1], as.character(pop(testset.gl))[1])
})

test_that("gl.write.csv works on SilicoDArT data", {
  td <- tempdir()
  gs <- testset.gs[1:5, 1:5]
  gs@other$loc.metrics <- testset.gs@other$loc.metrics[1:5, ]
  expect_no_error(
    capture.output(gl.write.csv(gs, outfile = "wcsv-gs.csv", outpath = td,
                                verbose = 0))
  )
  expect_true(file.exists(file.path(td, "wcsv-gs.csv")))
})

test_that("gl.write.csv returns NULL invisibly", {
  # [approved diff F1] no longer prints a bare NULL at the console
  td <- tempdir()
  v <- withVisible(
    gl.write.csv(testset.gl[1:5, 1:5], outfile = "wcsv-v.csv", outpath = td,
                 verbose = 0)
  )
  expect_null(v$value)
  expect_false(v$visible)
})

test_that("gl.write.csv falls back to tempdir for a non-existent outpath", {
  # [approved diff F2] gl.check.wd() now resolves outpath; a bad path falls
  # back to tempdir instead of throwing a connection error
  td <- tempdir()
  bad <- file.path(td, "wcsv-no-such-dir")
  expect_no_error(
    capture.output(gl.write.csv(testset.gl[1:5, 1:5], outfile = "wcsv-fb.csv",
                                outpath = bad, verbose = 0))
  )
  expect_false(file.exists(file.path(bad, "wcsv-fb.csv")))
  expect_true(file.exists(file.path(td, "wcsv-fb.csv")))
})

test_that("gl.write.csv leaves the input untouched", {
  td <- tempdir()
  x <- testset.gl
  capture.output(gl.write.csv(x, outfile = "wcsv-u.csv", outpath = td,
                              verbose = 0))
  expect_identical(x, dartR.data::testset.gl)
})

test_that("gl.write.csv message wording at verbose 2", {
  td <- tempdir()
  out <- capture.output(
    gl.write.csv(testset.gl[1:5, 1:5], outfile = "wcsv-vb.csv", outpath = td,
                 verbose = 2)
  )
  expect_true(any(grepl("Writing records to", out)))
  expect_true(any(grepl("Completed", out)))
})
