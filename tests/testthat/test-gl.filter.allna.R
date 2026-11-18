test_that("gl.filter.allna removes exactly the all-NA loci (SNP)", {
  x <- testset.gl
  m <- as.matrix(x)
  all_na <- apply(m, 2, function(col) all(is.na(col)))
  x2 <- gl.filter.allna(x, verbose = 0)
  expect_equal(nLoc(x2), nLoc(x) - sum(all_na))
  expect_setequal(locNames(x2), locNames(x)[!all_na])
  expect_equal(nrow(x2@other$loc.metrics), nLoc(x2))
  expect_equal(nInd(x2), nInd(x))
  expect_true(all(ploidy(x2) == 2))
  # history/flag updates only when something is removed (testset.gl in
  # dartR.data 1.2.2 has no all-NA loci, so this is a no-op call)
  if (sum(all_na) > 0) {
    expect_length(x2@other$history, length(x@other$history) + 1)
    expect_false(x2@other$loc.metrics.flags$CallRate)
  } else {
    expect_length(x2@other$history, length(x@other$history))
  }
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

# ---------------------------------------------------------------------------
# File-backed (FBM) objects (2026-11-10). All-missing loci and individuals
# are found from genotype counts read from the backing file in blocks
# (bigstatsr::big_counts), where as.matrix() read the whole matrix into
# memory and by.pop = TRUE wrote a copy of each population to disk first.
# The object returned must be the one the in-memory code returns.

fbm_allna_case <- function() {
  set.seed(3)
  m <- matrix(sample(0:2, 12 * 30, replace = TRUE), 12, 30,
              dimnames = list(sprintf("ind%02d", 1:12), sprintf("L%02d", 1:30)))
  m[sample.int(length(m), 40)] <- NA
  m[, 3] <- NA        # missing in every individual
  m[5:8, 7] <- NA     # missing in all of population B
  m[1:4, 14] <- NA    # missing in all of population A
  m[c(5, 12), ] <- NA # two individuals with no genotype at all
  x <- new("genlight", gen = m, ploidy = 2)
  pop(x) <- factor(rep(c("A", "B", "C"), each = 4))
  gl.compliance.check(x, verbose = 0)
}

expect_same_filtered <- function(b, a) {
  expect_identical(locNames(b), locNames(a))
  expect_identical(indNames(b), indNames(a))
  expect_equal(unname(as.matrix(b)), unname(as.matrix(a)))
  expect_identical(b@other$loc.metrics, a@other$loc.metrics)
  expect_identical(b@other$loc.metrics.flags, a@other$loc.metrics.flags)
  expect_identical(pop(b), pop(a))
  expect_false(is.null(b@fbm))
}

test_that("FBM-backed object: same loci and individuals removed as in memory", {
  skip_if_not_installed("bigstatsr")
  x <- fbm_allna_case()
  fb <- gl.gen2fbm(x, verbose = 0)
  a <- gl.filter.allna(x, verbose = 0)
  b <- gl.filter.allna(fb, verbose = 0)
  expect_equal(setdiff(locNames(x), locNames(b)), "L03")
  expect_setequal(setdiff(indNames(x), indNames(b)), c("ind05", "ind12"))
  expect_same_filtered(b, a)
  expect_length(b@other$history, length(fb@other$history) + 1)
})

test_that("FBM-backed object: by.pop = TRUE removes the same loci as in memory", {
  skip_if_not_installed("bigstatsr")
  x <- fbm_allna_case()
  # with individuals 5 and 12 present, B and C have loci missing by chance
  # too; the comparison is with the in-memory result, whatever it removes
  a <- gl.filter.allna(x, by.pop = TRUE, verbose = 0)
  b <- gl.filter.allna(gl.gen2fbm(x, verbose = 0), by.pop = TRUE, verbose = 0)
  expect_true(all(c("L03", "L07", "L14") %in% setdiff(locNames(x), locNames(b))))
  expect_same_filtered(b, a)

  a2 <- gl.filter.allna(testset.gl, by.pop = TRUE, verbose = 0)
  b2 <- gl.filter.allna(gl.gen2fbm(testset.gl, verbose = 0), by.pop = TRUE,
                        verbose = 0)
  expect_lt(nLoc(b2), nLoc(testset.gl))
  expect_same_filtered(b2, a2)
})

test_that("FBM-backed object: recalc = TRUE gives the same metrics as in memory", {
  skip_if_not_installed("bigstatsr")
  x <- fbm_allna_case()
  a <- gl.filter.allna(x, recalc = TRUE, verbose = 0)
  b <- gl.filter.allna(gl.gen2fbm(x, verbose = 0), recalc = TRUE, verbose = 0)
  expect_same_filtered(b, a)
})

test_that("FBM-backed object: nothing to remove returns the same data, no history entry", {
  skip_if_not_installed("bigstatsr")
  fb <- gl.gen2fbm(possums.gl, verbose = 0)
  b <- gl.filter.allna(fb, verbose = 0)
  expect_identical(locNames(b), locNames(fb))
  expect_identical(indNames(b), indNames(fb))
  expect_length(b@other$history, length(fb@other$history))
})
