# Characterization tests for gl2demerelate — snapshot of CURRENT behaviour
# (review baseline, function-review campaign). Bugs are captured as-is and
# annotated; do not treat every expectation here as intended behaviour.

test_that("gl2demerelate returns a Demerelate-shaped data frame for testset.gl", {
  df <- gl2demerelate(testset.gl, verbose = 0)

  expect_s3_class(df, "data.frame")
  expect_equal(dim(df), c(nInd(testset.gl), 2 + 2 * nLoc(testset.gl)))
  expect_equal(names(df)[1:2], c("Sample-ID", "Population"))
  # locus names are tidied (- and | to _, / removed) and suffixed _1/_2
  expect_equal(names(df)[3:4],
               c("100049687_12_CT_1", "100049687_12_CT_2"))
  expect_type(df[[1]], "character")
  expect_s3_class(df[[2]], "factor")
  # first individual, first locus: score 2 -> alleles (2,2); observed 2026-09-05
  expect_equal(unname(unlist(df[1, 3:4])), c(2, 2))
})

test_that("gl2demerelate allele coding: 0->(1,1), 1->(1,2), 2->(2,2), NA stays NA", {
  df <- gl2demerelate(testset.gl, verbose = 0)
  m <- as.matrix(testset.gl)

  tidy <- function(nm) {
    nm <- gsub("-", "_", nm, fixed = TRUE)
    nm <- gsub("/", "", nm, fixed = TRUE)
    gsub("|", "_", nm, fixed = TRUE)
  }
  check_one <- function(score, a1, a2) {
    idx <- which(m == score, arr.ind = TRUE)[1, ]
    loc <- tidy(colnames(m)[idx[2]])
    expect_equal(df[idx[1], paste0(loc, "_1")], a1)
    expect_equal(df[idx[1], paste0(loc, "_2")], a2)
  }
  check_one(0, 1, 1)
  check_one(1, 1, 2)
  check_one(2, 2, 2)

  naidx <- which(is.na(m), arr.ind = TRUE)[1, ]
  locna <- tidy(colnames(m)[naidx[2]])
  expect_true(is.na(df[naidx[1], paste0(locna, "_1")]))
  expect_true(is.na(df[naidx[1], paste0(locna, "_2")]))
})

test_that("gl2demerelate keeps the _1/_2 columns of each locus adjacent", {
  df <- gl2demerelate(testset.gl[1:10, 1:30], verbose = 0)
  loccols <- names(df)[-(1:2)]
  stems <- sub("_[12]$", "", loccols)
  expect_true(all(stems[seq(1, length(stems), 2)] ==
                  stems[seq(2, length(stems), 2)]))
})

test_that("gl2demerelate is silent at verbose = 0", {
  out <- capture.output(df <- gl2demerelate(testset.gl[1:5, 1:5], verbose = 0))
  expect_length(out, 0)
})

test_that("gl2demerelate rejects SilicoDArT data", {
  # [approved F1] accept = "SNP" now stops presence/absence scores from
  # being recoded as fake diploid genotypes (0 -> 1/1, 1 -> 1/2), which is
  # meaningless for Demerelate.
  expect_error(gl2demerelate(testset.gs[1:5, 1:10], verbose = 0),
               "SNP")
})
