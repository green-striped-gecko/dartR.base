# Chain-context characterization tests for utils.collapse.matrix
# Companion to test-utils.collapse.matrix.R (infrastructure-wave
# baseline, PR #325). Added during the population-distance chain review:
# pins the collapse statistic, label safety under non-alphabetical pop
# levels, NA propagation, and edge behaviour, at the post-#325 state on
# this branch. Assertions tagged [approved Fn] were flipped when the
# matching approved v2 finding was applied (PR #325 addendum).

fixture6 <- function() {
  gl6 <- testset.gl[1:6, 1:30]
  g <- new("genlight", gen = as.matrix(gl6), ploidy = 2)
  indNames(g) <- paste0("i", 1:6)
  locNames(g) <- locNames(gl6)
  # deliberately non-alphabetical level order
  pop(g) <- factor(c("A", "A", "A", "B", "B", "C"),
                   levels = c("B", "A", "C"))
  gl.compliance.check(g, verbose = 0)
}

fixtureM <- function() {
  M <- matrix(0, 6, 6,
              dimnames = list(paste0("i", 1:6), paste0("i", 1:6)))
  M[lower.tri(M)] <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
  M <- M + t(M)
  M["i1", "i2"] <- M["i2", "i1"] <- NA  # one NA pair within pop A
  M
}

test_that("between-pop cells are the mean over all cross-pop pairs", {
  g <- fixture6()
  M <- fixtureM()
  R <- as.matrix(utils.collapse.matrix(D = as.dist(M), x = g,
                                       verbose = 0))
  expect_equal(R["A", "B"],
               mean(M[c("i1", "i2", "i3"), c("i4", "i5")], na.rm = TRUE))
  expect_equal(R["A", "C"],
               mean(M[c("i1", "i2", "i3"), "i6"], na.rm = TRUE))
  expect_equal(R["B", "C"], mean(M[c("i4", "i5"), "i6"], na.rm = TRUE))
})

test_that("labels follow popNames order and are label-safe", {
  # non-alphabetical levels: name-based submatrix indexing keeps values
  # on the right labels (contrast with gl.dist.pop's dcast path)
  g <- fixture6()
  R <- as.matrix(utils.collapse.matrix(D = as.dist(fixtureM()), x = g,
                                       verbose = 0))
  expect_identical(rownames(R), c("B", "A", "C"))
})

test_that("matrix in, matrix out; diagonal over distinct pairs only", {
  # post-#325 state: within-pop mean excludes self-distances; a
  # single-individual population collapses to 0
  g <- fixture6()
  M <- fixtureM()
  R <- utils.collapse.matrix(D = M, x = g, verbose = 0)
  expect_true(is.matrix(R))
  expect_equal(R["A", "A"],
               mean(c(M["i1", "i3"], M["i2", "i3"]), na.rm = TRUE))
  expect_equal(R["B", "B"], M["i4", "i5"])
  expect_equal(R["C", "C"], 0)  # singleton population
})

test_that("dist in, dist out; the diagonal is dropped", {
  g <- fixture6()
  R <- utils.collapse.matrix(D = as.dist(fixtureM()), x = g, verbose = 0)
  expect_s3_class(R, "dist")
  expect_true(all(diag(as.matrix(R)) == 0))
})

test_that("an all-NA between-pop block propagates as NaN", {
  g <- fixture6()
  M <- fixtureM()
  M[c("i1", "i2", "i3"), "i6"] <- NA
  M["i6", c("i1", "i2", "i3")] <- NA
  R <- utils.collapse.matrix(D = M, x = g, verbose = 0)
  expect_true(is.nan(R["A", "C"]))
})

test_that("a D matrix missing an individual fails with a clear message", {
  # [approved F1] the name guard is now two-directional; previously a
  # D computed on a subset of individuals passed the one-way guard and
  # died later with a bare "subscript out of bounds"
  g <- fixture6()
  M <- fixtureM()[-6, -6]
  expect_error(utils.collapse.matrix(D = M, x = g, verbose = 0),
               "missing from the matrix")
})

test_that("an unnamed matrix fails with a clear message", {
  # [approved F2, folded into change 1] previously died at the %in%
  # guard with "no 'dimnames' attribute for array"
  g <- fixture6()
  M <- unname(fixtureM())
  expect_error(utils.collapse.matrix(D = M, x = g, verbose = 0),
               "must have row and column names")
})

test_that("verbose=0 is silent", {
  g <- fixture6()
  o <- capture.output(invisible(
    utils.collapse.matrix(D = as.dist(fixtureM()), x = g, verbose = 0)))
  expect_equal(length(o), 0)
})
