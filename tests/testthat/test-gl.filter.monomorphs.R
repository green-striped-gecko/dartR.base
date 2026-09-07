test_that("gl.filter.monomorphs removes exactly the monomorphic + all-NA loci (SNP)", {
  x <- testset.gl
  m <- as.matrix(x)
  is_mono <- apply(m, 2, function(col) {
    v <- col[!is.na(col)]
    length(v) == 0 || all(v == 0) || all(v == 2)
  })
  x2 <- gl.filter.monomorphs(x, verbose = 0)
  expect_equal(nLoc(x2), nLoc(x) - sum(is_mono))
  expect_setequal(locNames(x2), locNames(x)[!is_mono])
  expect_equal(nrow(x2@other$loc.metrics), nLoc(x2))
  expect_equal(nInd(x2), nInd(x))
  expect_true(all(ploidy(x2) == 2))
  expect_true(x2@other$loc.metrics.flags$monomorphs)
  expect_length(x2@other$history, length(x@other$history) + 1)
})

test_that("gl.filter.monomorphs removes exactly the monomorphic + all-NA loci (SilicoDArT)", {
  x <- testset.gs
  m <- as.matrix(x)
  is_mono <- apply(m, 2, function(col) {
    v <- col[!is.na(col)]
    length(v) == 0 || all(v == 0) || all(v == 1)
  })
  x2 <- gl.filter.monomorphs(x, verbose = 0)
  expect_equal(nLoc(x2), nLoc(x) - sum(is_mono))
  expect_setequal(locNames(x2), locNames(x)[!is_mono])
  expect_true(all(ploidy(x2) == 1))
})

test_that("all-heterozygous loci (uniform dosage 1) are correctly retained as polymorphic", {
  m2 <- matrix(c(1, 0, 1, 1, 2, 1, 1, 1, 1, 1, 0, 1), nrow = 3, byrow = TRUE,
               dimnames = list(paste0("I", 1:3), paste0("L", 1:4)))
  # L1: 1,1,1 all het -> polymorphic (retain); L2: 0,2,1 -> polymorphic;
  # L3: 1,1,0 -> polymorphic; L4: 1,1,1 all het -> retain
  gl_small <- new("genlight", gen = m2, ploidy = 2)
  pop(gl_small) <- factor(rep("p1", 3))
  gl_small <- gl.compliance.check(gl_small, verbose = 0)
  x2 <- gl.filter.monomorphs(gl_small, verbose = 0)
  expect_equal(nLoc(x2), 4)
})

test_that("a run with no monomorphs present returns the object with flag TRUE", {
  x <- gl.filter.monomorphs(testset.gl, verbose = 0)
  x2 <- gl.filter.monomorphs(x, verbose = 0)
  expect_equal(nLoc(x2), nLoc(x))
  expect_true(x2@other$loc.metrics.flags$monomorphs)
})

test_that("verbose = 0 is fully silent, including on an unassigned call (F1 fix)", {
  x <- testset.gl
  out0 <- capture.output(gl.filter.monomorphs(x, verbose = 0))
  expect_length(out0, 0)
})

test_that("history does not carry the internal gl.drop.loc entry (F2 fix)", {
  x <- testset.gl
  x2 <- gl.filter.monomorphs(x, verbose = 0)
  new_entries <- x2@other$history[(length(x@other$history) + 1):length(x2@other$history)]
  expect_length(new_entries, 1)
  expect_true(grepl("gl.filter.monomorphs", deparse(new_entries[[1]])[1]))
  expect_false(any(grepl("gl.drop.loc", vapply(new_entries, function(h) deparse(h)[1], ""))))
})
