# Characterization tests for gl2hiphop — snapshot of CURRENT behaviour
# (review baseline, bugs included). Captured 2026-09-05 at ed99203.

test_that("gl2hiphop recodes SNP genotypes to hiphop coding (0 hom, 1 other hom, 2 het)", {
  h <- gl2hiphop(testset.gl, verbose = 0)

  expect_s3_class(h, "data.frame")
  expect_equal(dim(h), c(274, 755))            # all-NA loci are retained

  # individuals as rownames, loci as colnames — the layout hiphop keys on
  expect_identical(rownames(h), indNames(testset.gl))
  expect_identical(colnames(h), locNames(testset.gl))

  m <- as.matrix(testset.gl)
  hm <- as.matrix(h)
  # dartR 0 (hom ref) -> 0; dartR 2 (hom alt) -> 1; dartR 1 (het) -> 2
  expect_equal(unique(hm[m == 0 & !is.na(m)]), 0)
  expect_equal(unique(hm[m == 2 & !is.na(m)]), 1)
  expect_equal(unique(hm[m == 1 & !is.na(m)]), 2)
  # NA survives as NA
  expect_identical(unname(is.na(m)), unname(is.na(hm)))
  expect_true(all(hm %in% c(0, 1, 2) | is.na(hm)))
})

test_that("gl2hiphop rejects SilicoDArT input", {
  # Expectation updated per approved finding F1 (DAT7): accept = 'SNP'
  # now rejects presence/absence data, which was previously recoded to
  # 0/2 pseudo-genotypes ('heterozygote' in hiphop coding) without any
  # warning.
  expect_error(gl2hiphop(testset.gs, verbose = 0), "SilicoDArT")
})

test_that("gl2hiphop is silent at verbose = 0", {
  out <- capture.output(h <- gl2hiphop(testset.gl, verbose = 0))
  expect_length(out, 0)
})
