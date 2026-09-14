# Characterization tests for gl2bayesAss — snapshot of CURRENT behaviour
# (review baseline, function-review campaign). Bugs are captured as-is and
# annotated; do not treat every expectation here as intended behaviour.

test_that("gl2bayesAss writes BA3 rows and returns them as a data.table (SNP)", {
  od <- file.path(tempdir(), "ba3_base")
  dir.create(od, showWarnings = FALSE)

  res <- gl2bayesAss(testset2.gl, outpath = od, verbose = 0)

  expect_s3_class(res, "data.table")
  expect_equal(dim(res), c(nInd(testset2.gl) * nLoc(testset2.gl), 5))
  expect_equal(names(res), c("rn", "Pop", "Locus", "All1", "All2"))

  f <- file.path(od, "gl.BayesAss.txt")
  expect_true(file.exists(f))
  ln <- readLines(f)
  expect_length(ln, 206870)  # 274 ind x 755 loci
  # sorted by Pop then individual; first line observed 2026-09-05
  expect_equal(ln[1], "AA033575 EmmacBrisWive 100049687-12-C/T 2 2")
})

test_that("gl2bayesAss genotype-to-allele recode: 0->(1,1) 1->(1,2) 2->(2,2) NA->(0,0)", {
  od <- file.path(tempdir(), "ba3_map")
  dir.create(od, showWarnings = FALSE)
  x <- testset2.gl[1:10, 1:20]
  res <- gl2bayesAss(x, outfile = "map.txt", outpath = od, verbose = 0)

  m <- as.matrix(x)
  key <- res[, paste(All1, All2)]
  geno <- m[cbind(match(res$rn, rownames(m)), match(res$Locus, colnames(m)))]
  expect_true(all(key[is.na(geno)] == "0 0"))
  expect_true(all(key[!is.na(geno) & geno == 0] == "1 1"))
  expect_true(all(key[!is.na(geno) & geno == 1] == "1 2"))
  expect_true(all(key[!is.na(geno) & geno == 2] == "2 2"))
})

test_that("gl2bayesAss is silent at verbose = 0", {
  od <- file.path(tempdir(), "ba3_sil")
  dir.create(od, showWarnings = FALSE)
  out <- capture.output(
    invisible(gl2bayesAss(testset2.gl[1:5, 1:5], outfile = "sil.txt",
                          outpath = od, verbose = 0)))
  expect_length(out, 0)
})

test_that("gl2bayesAss rejects ploidy != 2 with the house fatal-error idiom", {
  # [approved F3] plain stop() replaced by stop(error("Fatal Error: ..."))
  expect_error(gl2bayesAss(testset2.gl, ploidy = 4, verbose = 0),
               "only caters for ploidy = 2")
})

test_that("gl2bayesAss rejects SilicoDArT data", {
  # [approved F1] accept = "SNP" now stops presence/absence data from being
  # pushed through the diploid recode table.
  expect_error(gl2bayesAss(testset2.gs[1:5, 1:5], outfile = "gs.txt",
                           outpath = tempdir(), verbose = 0),
               "SNP")
})
