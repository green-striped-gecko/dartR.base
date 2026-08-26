test_that("gl.drop.ind removes the nominated individuals (SNP)", {
  x <- testset.gl
  to_drop <- c("AA019073", "AA004859")
  n_before <- nInd(x)

  x2 <- gl.drop.ind(x, ind.list = to_drop, verbose = 0)

  expect_false(any(to_drop %in% indNames(x2)))
  expect_equal(nInd(x2), n_before - length(to_drop))
  expect_equal(nLoc(x2), nLoc(x))
  expect_true(all(ploidy(x2) == 2))
  expect_length(x2@other$history, length(x@other$history) + 1)
})

test_that("gl.drop.ind works on SilicoDArT data", {
  x <- testset.gs
  to_drop <- indNames(x)[1:2]
  x2 <- gl.drop.ind(x, ind.list = to_drop, verbose = 0)
  expect_false(any(to_drop %in% indNames(x2)))
  expect_equal(nInd(x2), nInd(x) - 2)
})

test_that("gl.drop.ind mono.rm=TRUE + recalc=TRUE removes monomorphs and recalculates", {
  x <- testset.gl
  x2 <- gl.drop.ind(x, ind.list = c("AA019073", "AA004859"),
                     mono.rm = TRUE, recalc = TRUE, verbose = 0)
  m <- as.matrix(x2)
  is_mono <- apply(m, 2, function(col) length(unique(col[!is.na(col)])) <= 1)
  expect_false(any(is_mono))
  flg <- x2@other$loc.metrics.flags
  expect_true(all(unlist(flg[, sapply(flg, is.logical)])))
})

test_that("gl.drop.ind ignores a nonexistent individual name with a warning, not fatal", {
  x <- testset.gl
  real_ind <- indNames(x)[1]
  expect_error(
    x2 <- gl.drop.ind(x, ind.list = c(real_ind, "NOT_A_REAL_IND"), verbose = 0),
    NA
  )
  expect_false(real_ind %in% indNames(x2))
})

test_that("gl.drop.ind errors fatally when no listed individual exists", {
  x <- testset.gl
  expect_error(
    gl.drop.ind(x, ind.list = "NOT_A_REAL_IND", verbose = 0),
    "no individuals to drop"
  )
})

test_that("recalc=FALSE resets locus-metric flags to FALSE regardless of verbose (F1 fix)", {
  x <- testset.gl
  x0 <- gl.drop.ind(x, ind.list = c("AA019073", "AA004859"), recalc = FALSE, verbose = 0)
  flg0 <- x0@other$loc.metrics.flags
  expect_true(all(!unlist(flg0[, sapply(flg0, is.logical)])))

  x2 <- gl.drop.ind(x, ind.list = c("AA019073", "AA004859"), recalc = FALSE, verbose = 2)
  flg2 <- x2@other$loc.metrics.flags
  expect_true(all(!unlist(flg2[, sapply(flg2, is.logical)])))
})
