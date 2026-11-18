test_that("gl.keep.ind retains only the nominated individuals (SNP)", {
  x <- testset.gl
  to_keep <- c("AA019073", "AA004859")

  x2 <- gl.keep.ind(x, ind.list = to_keep, verbose = 0)

  expect_setequal(indNames(x2), to_keep)
  expect_equal(nInd(x2), length(to_keep))
  expect_equal(nLoc(x2), nLoc(x))
  expect_true(all(ploidy(x2) == 2))
  expect_length(x2@other$history, length(x@other$history) + 1)
})

test_that("gl.keep.ind works on SilicoDArT data", {
  x <- testset.gs
  to_keep <- indNames(x)[1:2]
  x2 <- gl.keep.ind(x, ind.list = to_keep, verbose = 0)
  expect_setequal(indNames(x2), to_keep)
})

test_that("gl.keep.ind mono.rm=TRUE + recalc=TRUE removes true monomorphs and recalculates", {
  x <- testset.gl
  x2 <- gl.keep.ind(x, ind.list = indNames(x)[1:20],
                     mono.rm = TRUE, recalc = TRUE, verbose = 0)
  m <- as.matrix(x2)
  # A locus is monomorphic only if every non-missing individual has the SAME
  # dosage AND that dosage is 0 or 2 (uniform dosage=1 means every individual
  # is heterozygous -- both alleles present, i.e. maximally polymorphic, not
  # monomorphic; see gl.keep.pop review for the false-positive this caused).
  is_mono <- apply(m, 2, function(col) {
    v <- unique(col[!is.na(col)])
    length(v) == 1 && v[1] %in% c(0, 2)
  })
  expect_false(any(is_mono))
  flg <- x2@other$loc.metrics.flags
  expect_true(all(unlist(flg[, sapply(flg, is.logical)])))
})

test_that("gl.keep.ind ignores a nonexistent individual name with a warning, not fatal", {
  x <- testset.gl
  real_ind <- indNames(x)[1]
  expect_error(
    x2 <- gl.keep.ind(x, ind.list = c(real_ind, "NOT_A_REAL_IND"), verbose = 0),
    NA
  )
  expect_true(real_ind %in% indNames(x2))
})

test_that("gl.keep.ind errors fatally when no listed individual exists", {
  x <- testset.gl
  expect_error(
    gl.keep.ind(x, ind.list = "NOT_A_REAL_IND", verbose = 0),
    "no individuals listed to keep"
  )
})

test_that("recalc=FALSE resets locus-metric flags to FALSE regardless of verbose (F1 fix)", {
  x <- testset.gl
  x0 <- gl.keep.ind(x, ind.list = indNames(x)[1:20], recalc = FALSE, verbose = 0)
  flg0 <- x0@other$loc.metrics.flags
  expect_true(all(!unlist(flg0[, sapply(flg0, is.logical)])))

  x2 <- gl.keep.ind(x, ind.list = indNames(x)[1:20], recalc = FALSE, verbose = 2)
  flg2 <- x2@other$loc.metrics.flags
  expect_true(all(!unlist(flg2[, sapply(flg2, is.logical)])))
})
