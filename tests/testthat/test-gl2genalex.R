# Characterization tests for gl2genalex — snapshot of CURRENT behaviour (bugs
# included) captured pre-review at commit ed99203. Not a statement of correct
# behaviour; accepted diffs must map to approved review findings.

test_that("gl2genalex writes a GenAlEx codominant csv with 1-4 allele codes and 0 missing", {
  skip_if_not_installed("poppr")
  small <- gl.filter.allna(testset2.gl[1:6, 1:8], verbose = 0)
  # small is 6 individuals x 7 loci, 6 populations
  td <- file.path(tempdir(), "galex_char")
  dir.create(td, showWarnings = FALSE)
  out <- capture.output(res <- gl2genalex(small, outfile = "small.csv",
                                          outpath = td, overwrite = TRUE,
                                          verbose = 0))
  expect_null(res)
  f <- file.path(td, "small.csv")
  expect_true(file.exists(f))
  lns <- readLines(f)
  # 3 header rows + one row per individual
  expect_length(lns, 3 + nInd(small))
  # header row 1: n loci, n individuals, n pops, then pop sizes
  h1 <- strsplit(lns[1], ",")[[1]]
  expect_identical(h1[1:9], c("7", "6", "6", "1", "1", "1", "1", "1", "1"))
  # header row 2 carries pop names from column 4
  h2 <- strsplit(lns[2], ",")[[1]]
  expect_identical(h2[4:9], c("EmmacBurnBara", "EmmacCoopEulb", "EmmacMaclGeor",
                              "EmmacMDBForb", "EmmacMDBMaci", "EmmacMDBSanf"))
  # header row 3: Ind, Pop, then locus names with a blank second-allele column
  expect_true(startsWith(lns[3], "Ind,Pop,100049687-12-C/T,"))
  # first data row snapshot: alleles coded A=1, C=2, G=3, T=4, missing 0
  expect_identical(
    lns[4],
    paste0("\"AA010915\",\"EmmacMDBForb\",\"4\",\"4\",\"1\",\"1\",\"4\",\"4\",",
           "\"1\",\"1\",\"0\",\"0\",\"2\",\"2\",\"3\",\"4\"")
  )
})

test_that("gl2genalex is silent at verbose = 0", {
  skip_if_not_installed("poppr")
  # Updated for approved findings F3 (quiet = (verbose < 2) passed to
  # poppr::genind2genalex) and F5 (invisible NULL return).
  small <- gl.filter.allna(testset2.gl[1:6, 1:8], verbose = 0)
  td <- file.path(tempdir(), "galex_char")
  dir.create(td, showWarnings = FALSE)
  out <- capture.output(vis <- withVisible(
    gl2genalex(small, outfile = "v0.csv", outpath = td,
               overwrite = TRUE, verbose = 0)))
  expect_false(any(grepl("Extracting the table", out)))
  expect_length(out, 0)
  expect_null(vis$value)
  expect_false(vis$visible)
})

test_that("gl2genalex overwrite = FALSE leaves an existing file in place", {
  skip_if_not_installed("poppr")
  small <- gl.filter.allna(testset2.gl[1:6, 1:8], verbose = 0)
  td <- file.path(tempdir(), "galex_char")
  dir.create(td, showWarnings = FALSE)
  f <- file.path(td, "keep.csv")
  writeLines("sentinel", f)
  # current behaviour: poppr errors when the file exists and overwrite = FALSE
  expect_error(
    capture.output(gl2genalex(small, outfile = "keep.csv", outpath = td,
                              overwrite = FALSE, verbose = 0))
  )
  expect_identical(readLines(f), "sentinel")
})

test_that("gl2genalex rejects SilicoDArT input", {
  skip_if_not_installed("poppr")
  # Updated for approved finding F1 (DAT7): accept = "SNP" now rejects
  # presence/absence data instead of exporting pseudo-codominant genotypes.
  gs <- testset2.gs[1:4, 1:6]
  gs@other$loc.metrics <- testset2.gs@other$loc.metrics[1:6, ]
  td <- file.path(tempdir(), "galex_char")
  dir.create(td, showWarnings = FALSE)
  expect_error(
    capture.output(gl2genalex(gs, outfile = "gs.csv", outpath = td,
                              overwrite = TRUE, verbose = 0)),
    "inappropriate object passed to function"
  )
})
