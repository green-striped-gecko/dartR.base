test_that("gl.filter.allna removes exactly the all-NA loci (SNP, testset has 3)", {
  x <- testset.gl
  m <- as.matrix(x)
  all_na <- apply(m, 2, function(col) all(is.na(col)))
  x2 <- gl.filter.allna(x, verbose = 0)
  expect_equal(nLoc(x2), nLoc(x) - sum(all_na))
  expect_setequal(locNames(x2), locNames(x)[!all_na])
  expect_equal(nrow(x2@other$loc.metrics), nLoc(x2))
  expect_equal(nInd(x2), nInd(x))
  expect_true(all(ploidy(x2) == 2))
  expect_length(x2@other$history, length(x@other$history) + 1)
  expect_false(x2@other$loc.metrics.flags$CallRate)
})

test_that("gl.filter.allna works on SilicoDArT", {
  x <- testset.gs
  m <- as.matrix(x)
  all_na <- apply(m, 2, function(col) all(is.na(col)))
  x2 <- gl.filter.allna(x, verbose = 0)
  expect_equal(nLoc(x2), nLoc(x) - sum(all_na))
  expect_true(all(ploidy(x2) == 1))
})

test_that("gl.filter.allna by.pop=TRUE removes the union of per-population all-NA loci", {
  x <- testset.gl
  expected <- unique(unlist(lapply(levels(pop(x)), function(p) {
    sub <- as.matrix(x[pop(x) == p, ])
    locNames(x)[apply(sub, 2, function(col) all(is.na(col)))]
  })))
  x2 <- gl.filter.allna(x, by.pop = TRUE, verbose = 0)
  expect_equal(nLoc(x2), nLoc(x) - length(expected))
  expect_false(any(expected %in% locNames(x2)))
})

test_that("gl.filter.allna removes an all-NA individual correctly", {
  m2 <- matrix(c(0, 1, 2, 0, NA, NA, NA, NA, 2, 1, 0, 1), nrow = 3, byrow = TRUE,
               dimnames = list(c("IND_A", "IND_ALLNA", "IND_C"), paste0("L", 1:4)))
  gl_small <- new("genlight", gen = m2, ploidy = 2)
  pop(gl_small) <- factor(rep("p1", nInd(gl_small)))
  gl_small <- gl.compliance.check(gl_small, verbose = 0)
  x2 <- gl.filter.allna(gl_small, verbose = 0)
  expect_equal(nInd(x2), 2)
  expect_false("IND_ALLNA" %in% indNames(x2))
  expect_equal(nLoc(x2), 4)
})

test_that("verbose = 0 is fully silent, including on an unassigned call (F2 fix)", {
  x <- testset.gl
  out0 <- capture.output(gl.filter.allna(x, verbose = 0))
  expect_length(out0, 0)
})

test_that("individuals-only removal invalidates flags and records history (F1 fix)", {
  m2 <- matrix(c(0, 1, 2, 0, NA, NA, NA, NA, 2, 1, 0, 1), nrow = 3, byrow = TRUE,
               dimnames = list(c("IND_A", "IND_ALLNA", "IND_C"), paste0("L", 1:4)))
  gl_small <- new("genlight", gen = m2, ploidy = 2)
  pop(gl_small) <- factor(rep("p1", nInd(gl_small)))
  gl_small <- gl.compliance.check(gl_small, verbose = 0)
  gl_small@other$loc.metrics.flags$CallRate <- TRUE
  h0 <- length(gl_small@other$history)
  x2 <- gl.filter.allna(gl_small, verbose = 0)
  expect_false("IND_ALLNA" %in% indNames(x2))
  expect_false(x2@other$loc.metrics.flags$CallRate)
  expect_equal(length(x2@other$history) - h0, 1)
})

test_that("all-NA individuals listing has a count and no NULLs (F5 fix)", {
  m2 <- matrix(c(0, 1, 2, 0, NA, NA, NA, NA, 2, 1, 0, 1), nrow = 3, byrow = TRUE,
               dimnames = list(c("IND_A", "IND_ALLNA", "IND_C"), paste0("L", 1:4)))
  gl_small <- new("genlight", gen = m2, ploidy = 2)
  pop(gl_small) <- factor(rep("p1", nInd(gl_small)))
  gl_small <- gl.compliance.check(gl_small, verbose = 0)
  outc <- capture.output(gl.filter.allna(gl_small, verbose = 3))
  ind_line <- grep("individuals that are missing", outc, value = TRUE)
  expect_length(ind_line, 1)
  expect_true(grepl("IND_ALLNA", ind_line))
  expect_false(any(grepl("NULL", outc)))
})

test_that("by.pop records one history entry, no internal gl.drop.loc leak (F4 fix)", {
  x <- testset.gl
  h1 <- length(x@other$history)
  xp <- gl.filter.allna(x, by.pop = TRUE, verbose = 0)
  expect_equal(length(xp@other$history) - h1, 1)
  new_entry <- xp@other$history[[length(xp@other$history)]]
  expect_true(grepl("gl.filter.allna", deparse(new_entry)[1]))
})

test_that("non-genlight input fails fast with the standard message (F3 fix)", {
  expect_error(
    capture.output(gl.filter.allna(data.frame(a = 1), verbose = 0)),
    "inappropriate object passed"
  )
})
