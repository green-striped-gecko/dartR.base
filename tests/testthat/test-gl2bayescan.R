# Characterization tests for gl2bayescan — snapshot of CURRENT behaviour
# (review baseline, function-review campaign). Bugs are captured as-is and
# annotated; do not treat every expectation here as intended behaviour.

test_that("gl2bayescan writes the BayeScan GESTE layout for testset.gl", {
  od <- file.path(tempdir(), "bsc_base")
  dir.create(od, showWarnings = FALSE)

  ret <- gl2bayescan(testset.gl, outfile = "bayescan.txt", outpath = od,
                     verbose = 0)
  expect_null(ret)

  f <- file.path(od, "bayescan.txt")
  expect_true(file.exists(f))
  bl <- readLines(f)
  expect_length(bl, 23471)
  expect_equal(trimws(bl[1]), "[loci]=755")
  expect_equal(trimws(bl[3]), "[populations]=31")
  expect_equal(length(grep("^\\[pop\\]", bl)), nPop(testset.gl))
})

test_that("gl2bayescan allele counts match manual counts for pop 1, locus 1", {
  od <- file.path(tempdir(), "bsc_cnt")
  dir.create(od, showWarnings = FALSE)
  invisible(gl2bayescan(testset.gl, outfile = "cnt.txt", outpath = od,
                        verbose = 0))
  bl <- readLines(file.path(od, "cnt.txt"))

  p1 <- as.character(unique(pop(testset.gl)))[1]
  msub <- as.matrix(testset.gl[pop(testset.gl) == p1, ])
  nobs1 <- sum(!is.na(msub[, 1]))
  alt1 <- sum(msub[, 1], na.rm = TRUE)

  first_data <- trimws(bl[grep("^\\[pop\\]", bl)[1] + 1])
  expect_equal(first_data,
               paste(1, 2 * nobs1, 2, alt1, 2 * nobs1 - alt1))
  # observed 2026-09-05: "1 22 2 22 0"
  expect_equal(first_data, "1 22 2 22 0")
})

test_that("gl2bayescan writes zero-sample rows for all-missing pop x locus cells (captured as-is)", {
  # Current behaviour: a locus with no genotyped individuals in a population
  # is written as "<j> 0 2 0 0". BayeScan's handling of zero gene copies is
  # not verified; captured so a fix (dropping or warning) is visible.
  od <- file.path(tempdir(), "bsc_zero")
  dir.create(od, showWarnings = FALSE)
  invisible(gl2bayescan(testset.gl, outfile = "zero.txt", outpath = od,
                        verbose = 0))
  bl <- readLines(file.path(od, "zero.txt"))
  datl <- grep("^[0-9]+ ", trimws(bl), value = TRUE)
  genes <- as.numeric(sapply(strsplit(datl, " +"), `[`, 2))
  expect_length(datl, nLoc(testset.gl) * nPop(testset.gl))
  expect_equal(sum(genes == 0), 815)

  # [approved F2] the zero-sample count is now reported at verbose >= 1
  wout <- capture.output(
    invisible(gl2bayescan(testset.gl, outfile = "zero2.txt", outpath = od,
                          verbose = 1)))
  expect_true(any(grepl("815 population x locus", wout)))
})

test_that("gl2bayescan is silent at verbose = 0 and returns NULL invisibly", {
  od <- file.path(tempdir(), "bsc_sil")
  dir.create(od, showWarnings = FALSE)
  out <- capture.output(
    v <- withVisible(gl2bayescan(testset.gl[1:5, 1:10], outfile = "sil.txt",
                                 outpath = od, verbose = 0)))
  expect_length(out, 0)
  # [approved F4] return(invisible(NULL)) -- an unassigned call no longer
  # prints NULL
  expect_false(v$visible)
})

test_that("gl2bayescan rejects SilicoDArT data", {
  # [approved F1] accept = "SNP" now stops presence/absence tag counts from
  # being written into a codominant (diploid gene-copy) BayeScan file.
  expect_error(gl2bayescan(testset.gs[1:5, 1:10], outfile = "gs.txt",
                           outpath = tempdir(), verbose = 0),
               "SNP")
})
