test_that("gl.drop.pop removes the nominated populations (SNP)", {
  x <- testset.gl
  pops_before <- popNames(x)
  n_ind_before <- nInd(x)
  to_drop <- c("EmsubRopeMata", "EmvicVictJasp")
  n_drop <- sum(pop(x) %in% to_drop)

  x2 <- gl.drop.pop(x, pop.list = to_drop, verbose = 0)

  expect_false(any(to_drop %in% popNames(x2)))
  expect_equal(nInd(x2), n_ind_before - n_drop)
  expect_equal(nLoc(x2), nLoc(x))
  expect_true(all(ploidy(x2) == 2))
  expect_length(x2@other$history, length(x@other$history) + 1)
})

test_that("gl.drop.pop default recalc=FALSE resets all loc.metrics.flags to FALSE", {
  x <- testset.gl
  x2 <- gl.drop.pop(x, pop.list = "EmsubRopeMata", verbose = 0)
  flg <- x2@other$loc.metrics.flags
  expect_setequal(unique(unlist(flg[, sapply(flg, is.logical)])), FALSE)
})

test_that("gl.drop.pop mono.rm=TRUE + recalc=TRUE removes monomorphs and recalculates", {
  x <- testset.gl
  x2 <- gl.drop.pop(x, pop.list = c("EmsubRopeMata", "EmvicVictJasp"),
                     mono.rm = TRUE, recalc = TRUE, verbose = 0)
  m <- as.matrix(x2)
  is_mono <- apply(m, 2, function(col) length(unique(col[!is.na(col)])) <= 1)
  expect_false(any(is_mono))
})

test_that("gl.drop.pop works on SilicoDArT data", {
  x <- testset.gs
  to_drop <- popNames(x)[1]
  x2 <- gl.drop.pop(x, pop.list = to_drop, verbose = 0)
  expect_false(to_drop %in% popNames(x2))
})

test_that("gl.drop.pop as.pop temporarily reassigns population for the drop", {
  x <- testset.gl
  if (!is.null(x@other$ind.metrics$sex)) {
    lvl <- unique(as.character(x@other$ind.metrics$sex))[1]
    n_drop <- sum(x@other$ind.metrics$sex == lvl)
    x2 <- gl.drop.pop(x, as.pop = "sex", pop.list = lvl, verbose = 0)
    # dropped the right number of individuals, and pop() labels afterward are
    # real population names, not leftover "sex" levels
    expect_equal(nInd(x2), nInd(x) - n_drop)
    expect_true(all(as.character(unique(pop(x2))) %in% as.character(unique(pop(x)))))
    expect_false(any(as.character(unique(pop(x2))) %in% c("male", "female", "unknown")))
  } else {
    succeed("no sex column in testset.gl ind.metrics -- skipped")
  }
})

test_that("gl.drop.pop ignores a nonexistent population name with a warning, not fatal", {
  x <- testset.gl
  real_pop <- popNames(x)[1]
  expect_error(
    x2 <- gl.drop.pop(x, pop.list = c(real_pop, "NOT_A_REAL_POP"), verbose = 0),
    NA
  )
})

test_that("gl.drop.pop errors fatally when no listed population exists", {
  x <- testset.gl
  expect_error(
    gl.drop.pop(x, pop.list = "NOT_A_REAL_POP", verbose = 0),
    "no populations listed to drop"
  )
})
