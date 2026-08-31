test_that("gl.keep.pop retains only the nominated populations (SNP)", {
  x <- testset.gl
  to_keep <- c("EmsubRopeMata", "EmvicVictJasp")
  n_keep <- sum(pop(x) %in% to_keep)

  x2 <- gl.keep.pop(x, pop.list = to_keep, verbose = 0)

  expect_setequal(popNames(x2), to_keep)
  expect_equal(nInd(x2), n_keep)
  expect_equal(nLoc(x2), nLoc(x))
  expect_true(all(ploidy(x2) == 2))
  expect_length(x2@other$history, length(x@other$history) + 1)
})

test_that("gl.keep.pop works on SilicoDArT data", {
  x <- testset.gs
  to_keep <- popNames(x)[1:2]
  x2 <- gl.keep.pop(x, pop.list = to_keep, verbose = 0)
  expect_true(all(popNames(x2) %in% to_keep))
})

test_that("gl.keep.pop mono.rm=TRUE + recalc=TRUE removes monomorphs and recalculates", {
  x <- testset.gl
  x2 <- gl.keep.pop(x, pop.list = c("EmsubRopeMata", "EmvicVictJasp"),
                     mono.rm = TRUE, recalc = TRUE, verbose = 0)
  m <- as.matrix(x2)
  is_mono <- apply(m, 2, function(col) {
    v <- unique(col[!is.na(col)])
    length(v) == 1 && v[1] %in% c(0, 2)
  })
  expect_false(any(is_mono))
  flg <- x2@other$loc.metrics.flags
  expect_true(all(unlist(flg[, sapply(flg, is.logical)])))
})

test_that("gl.keep.pop as.pop temporarily reassigns population for the keep", {
  x <- testset.gl
  if (!is.null(x@other$ind.metrics$sex)) {
    lvl <- unique(as.character(x@other$ind.metrics$sex))[1]
    n_keep <- sum(x@other$ind.metrics$sex == lvl)
    x2 <- gl.keep.pop(x, as.pop = "sex", pop.list = lvl, verbose = 0)
    expect_equal(nInd(x2), n_keep)
    expect_true(all(as.character(unique(pop(x2))) %in% as.character(unique(pop(x)))))
    expect_false(any(as.character(unique(pop(x2))) %in% c("male", "female", "unknown")))
  } else {
    succeed("no sex column in testset.gl ind.metrics -- skipped")
  }
})

test_that("gl.keep.pop ignores a nonexistent population name with a warning, not fatal", {
  x <- testset.gl
  real_pop <- popNames(x)[1]
  expect_error(
    x2 <- gl.keep.pop(x, pop.list = c(real_pop, "NOT_A_REAL_POP"), verbose = 0),
    NA
  )
  expect_true(real_pop %in% popNames(x2))
})

test_that("gl.keep.pop errors fatally when no listed population exists", {
  x <- testset.gl
  expect_error(
    gl.keep.pop(x, pop.list = "NOT_A_REAL_POP", verbose = 0),
    "no populations listed to keep"
  )
})

test_that("gl.keep.pop resets locus-metric flags to FALSE regardless of verbose (checking for the gl.drop.ind F1 pattern)", {
  x <- testset.gl
  x0 <- gl.keep.pop(x, pop.list = c("EmsubRopeMata", "EmvicVictJasp"), recalc = FALSE, verbose = 0)
  flg0 <- x0@other$loc.metrics.flags
  expect_true(all(!unlist(flg0[, sapply(flg0, is.logical)])))

  x2 <- gl.keep.pop(x, pop.list = c("EmsubRopeMata", "EmvicVictJasp"), recalc = FALSE, verbose = 2)
  flg2 <- x2@other$loc.metrics.flags
  expect_true(all(!unlist(flg2[, sapply(flg2, is.logical)])))
})
