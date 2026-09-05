# Characterization tests for gl2phylip
# Baseline snapshotted pre-review (dev at ed99203). Assertions capture
# CURRENT behaviour, defects included, so Phase C changes surface as diffs.

test_that("distance matrix and file structure are stable (bstrap = 1)", {
  td <- tempdir()
  d <- gl2phylip(testset.gl, outfile = "chphy1.txt", outpath = td, verbose = 0)
  expect_true(is.matrix(d))
  expect_equal(dim(d), c(31, 31))
  # names padded/truncated to the phylip 10-character width
  expect_equal(unique(nchar(rownames(d))), 10)
  expect_equal(sum(is.na(d)), 0)
  # F3: frequencies now computed with na.rm = TRUE, so distances change on
  # data with missing genotypes (values recomputed under the fix)
  expect_equal(unname(d[1, 2]), 7.7752)
  ln <- readLines(file.path(td, "chphy1.txt"))
  expect_length(ln, 32)
  expect_equal(ln[1], "    31 ")
  expect_match(ln[2], "^EmmacBrisW 0 7\\.7752 7\\.6009")
})

test_that("bstrap > 1 appends replicate matrices and returns the observed matrix", {
  td <- tempdir()
  d1 <- gl2phylip(testset.gl, outfile = "chphyA.txt", outpath = td, verbose = 0)
  set.seed(1)
  d3 <- gl2phylip(testset.gl, outfile = "chphy3.txt", outpath = td,
                  bstrap = 3, verbose = 0)
  ln <- readLines(file.path(td, "chphy3.txt"))
  expect_length(ln, 96)  # three 32-line matrices
  # F1: the caller now receives the observed-data matrix, not a replicate
  expect_true(isTRUE(all.equal(unname(d3), unname(d1))))
})

test_that("verbose = 0 is silent and the sink is released", {
  td <- tempdir()
  o <- capture.output(invisible(gl2phylip(testset.gl, outfile = "chphyv.txt",
                                          outpath = td, verbose = 0)))
  expect_length(o, 0)
  expect_equal(sink.number(), 0)
})

test_that("SilicoDArT input is rejected", {
  # F4: the diploid-dosage math is SNP-only; SilicoDArT now errors at the
  # datatype gate instead of converting with halved frequencies
  gs <- testset.gs[1:20, 1:100]
  gs@other$loc.metrics <- gs@other$loc.metrics[1:100, ]
  td <- tempdir()
  expect_error(capture.output(gl2phylip(gs, outfile = "chphygs.txt",
        outpath = td, verbose = 0)), "expecting SNP")
})
