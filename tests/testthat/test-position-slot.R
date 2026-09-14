# Characterization tests for the genome-only @position design
# (branch position-genome-only, dev at ddaed27). The @position/@chromosome
# slots are reserved for genome coordinates; the position of the SNP
# within the sequence tag lives in @other$loc.metrics$SnpPosition
# (0-based). DArT-read objects carry @position = NULL until genome
# coordinates are assigned explicitly.

test_that("the DArT reader leaves @position NULL and keeps SnpPosition", {
  csv <- system.file("extdata", "testset_SNPs_2Row.csv", package = "dartR.data")
  skip_if(csv == "", "dartR.data extdata not available")
  o <- capture.output(suppressWarnings(
    x <- gl.read.dart(csv, verbose = 0)))
  expect_null(x@position)
  expect_null(x@chromosome)
  sp <- x@other$loc.metrics$SnpPosition
  expect_equal(length(sp), nLoc(x))
  expect_true(min(as.numeric(as.character(sp)), na.rm = TRUE) >= 0)
})

test_that("gl.compliance.check clears a provably stale tag-position copy", {
  s <- testset.gl[1:6, 1:12]
  expect_identical(as.integer(s@position),
                   as.integer(s@other$loc.metrics$SnpPosition))
  o <- capture.output(s2 <- gl.compliance.check(s, verbose = 2))
  expect_null(s2@position)
  expect_true(any(grepl("stale copy of SnpPosition", o)))
  expect_equal(length(s2@other$loc.metrics$SnpPosition), nLoc(s2))
})

test_that("gl.compliance.check preserves genuine genome coordinates", {
  s <- testset.gl[1:6, 1:12]
  s@position <- as.integer(seq(1e6, 2e6, length.out = nLoc(s)))
  s@chromosome <- factor(rep("chr1", nLoc(s)))
  o <- capture.output(s2 <- gl.compliance.check(s, verbose = 0))
  expect_identical(s2@position, s@position)
  expect_identical(s2@chromosome, s@chromosome)
})

test_that("gl.compliance.check no longer refills a NULL @position", {
  s <- testset.gl[1:6, 1:12]
  s@position <- NULL
  o <- capture.output(s2 <- gl.compliance.check(s, verbose = 2))
  expect_null(s2@position)
  expect_false(any(grepl("Assigning SNP position", o)))
})

test_that("gl2bpp runs on a NULL-@position object via SnpPosition", {
  sub <- testset.gl[1:6, 1:12]
  o0 <- capture.output(sub <- gl.compliance.check(sub, verbose = 0))
  expect_null(sub@position)
  td <- tempdir()
  o <- capture.output(gl2bpp(sub, outfile = "pos_bpp.txt", outpath = td,
                             verbose = 0))
  expect_true(file.exists(file.path(td, "pos_bpp.txt")))
  ln <- readLines(file.path(td, "pos_bpp.txt"))
  expect_gt(length(ln), nInd(sub))
})

test_that("gl2plink's dummy fallback engages on a NULL-@position object", {
  sub <- testset.gl[1:6, 1:12]
  o0 <- capture.output(sub <- gl.compliance.check(sub, verbose = 0))
  td <- file.path(tempdir(), "posplink")
  dir.create(td, showWarnings = FALSE)
  o <- capture.output(suppressWarnings(
    gl2plink(sub, outfile = "pp.plink", outpath = td, verbose = 0)))
  mp <- list.files(td, pattern = "map$", full.names = TRUE)
  expect_true(length(mp) >= 1)
})

test_that("gl2hapmap uses genome coordinates when present, zeros when absent", {
  sub <- testset.gl[1:6, 1:12]
  o0 <- capture.output(sub <- gl.compliance.check(sub, verbose = 0))
  td <- file.path(tempdir(), "poshap")
  dir.create(td, showWarnings = FALSE)
  # absent -> zeroed, and no max() warning on the NULL slot
  expect_warning(
    o1 <- capture.output(h1 <- gl2hapmap(sub, outfile = "h1.hmp.txt",
                                         outpath = td, verbose = 0)),
    NA)
  # present -> passed through
  g <- sub
  g@position <- as.integer(seq(1e6, 2e6, length.out = nLoc(g)))
  g@chromosome <- factor(rep("1", nLoc(g)))
  o2 <- capture.output(h2 <- gl2hapmap(g, outfile = "h2.hmp.txt",
                                       outpath = td, verbose = 0))
  f2 <- list.files(td, pattern = "h2", full.names = TRUE)
  if (length(f2)) {
    hm <- readLines(f2[1])
    expect_true(any(grepl("1000000", hm)))
  }
})
