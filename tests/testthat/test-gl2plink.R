# Characterization tests for gl2plink
# Baseline snapshotted pre-review (dev at ed99203). Assertions capture
# CURRENT behaviour, defects included, so Phase C changes surface as diffs.
# bed.files = TRUE is not exercised (needs an external PLINK binary).

make_platy <- function() {
  test <- platypus.gl
  test$position <- test$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
  test$chromosome <- as.factor("1")
  test
}

test_that(".map and .ped structure with assigned position/chromosome", {
  td <- tempdir()
  t1 <- make_platy()
  o <- capture.output(rv <- gl2plink(t1, outfile = "chplk1", outpath = td,
                                     verbose = 0))
  expect_null(rv)
  map <- readLines(file.path(td, "chplk1.map"))
  expect_length(map, 1000)
  expect_equal(map[1], "1 45055704-36-G/A 0 2438118")
  ped <- readLines(file.path(td, "chplk1.ped"))
  expect_length(ped, 81)
  f <- strsplit(ped[1], " ")[[1]]
  expect_length(f, 2006)  # 6 leading columns + 2 alleles x 1000 loci
  # F2: the default scalar sex.code 'unknown' is now recoded to PLINK's 0
  expect_equal(f[1:6], c("TENTERFIELD", "T27", "0", "0", "0", "0"))
  # genotype coding: 0 -> hom ref alleles, NA -> 0 0
  expect_equal(f[7:8], c("G", "G"))
  expect_equal(f[11:12], c("0", "0"))  # locus 3 is NA for individual 1
})

test_that("sex.code is recoded for scalar and vector forms", {
  td <- tempdir()
  # F2: scalar sex.code is now recoded (recycling handles the scalar)
  capture.output(gl2plink(make_platy(), outfile = "chplkF", outpath = td,
                          sex.code = "F", verbose = 0))
  ped <- readLines(file.path(td, "chplkF.ped"))
  expect_equal(strsplit(ped[1], " ")[[1]][5], "2")
  capture.output(gl2plink(make_platy(), outfile = "chplkV", outpath = td,
                          sex.code = rep("F", nInd(platypus.gl)), verbose = 0))
  pedv <- readLines(file.path(td, "chplkV.ped"))
  expect_equal(strsplit(pedv[1], " ")[[1]][5], "2")
})

test_that("NULL position/chromosome are coded as unmapped (0/0) with gated warnings", {
  td <- tempdir()
  # F5: the fallback warnings are gated; verbose = 0 is silent
  o <- capture.output(gl2plink(testset.gl, outfile = "chplkT", outpath = td,
                               verbose = 0))
  expect_length(o, 0)
  # F5: at verbose >= 1 the fallbacks warn
  o1 <- capture.output(gl2plink(testset.gl, outfile = "chplkT", outpath = td,
                                verbose = 1))
  expect_true(any(grepl("chromosome slot is empty", o1)))
  expect_true(any(grepl("position slot is empty", o1)))
  map <- readLines(file.path(td, "chplkT.map"))
  expect_length(map, 755)
  # F1: unmapped loci follow PLINK's convention, chromosome 0 / position 0,
  # instead of fabricated chromosome 1 and sequential positions
  expect_equal(map[1], "0 100049687-12-C/T 0 0")
  expect_equal(map[755], "0 SIM500-46-C/G 0 0")
})

test_that("SilicoDArT input is rejected at the datatype gate", {
  # F3: accept = "SNP" replaces the mid-write crash with a clear error
  gs <- testset.gs[1:20, 1:100]
  gs@other$loc.metrics <- gs@other$loc.metrics[1:100, ]
  expect_error(capture.output(gl2plink(gs, outfile = "chplkGS",
      outpath = tempdir(), verbose = 0)),
      "expecting SNP")
})
