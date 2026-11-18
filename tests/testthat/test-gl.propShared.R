# Characterization tests for gl.propShared
# Baseline snapshotted before review (dartR.base dev at ddaed27), re-run
# unchanged on dartR.spatial dev at 3bba9b7. Assertions tagged [fixed: n]
# were flipped by approved change n of
# function-review/reports/dartR.spatial/gl.propShared.md.
#
# The equivalence pins in the last block compare against gl.dist.ind, whose
# SNP engine (utils.dist.ind.snp) carries its applied PR #315 fixes on this
# branch. At ddaed27 as released, only method = "manhattan"/"czekanowski"
# matches; "simple" was reference-allele asymmetric until #315.

test_that("returned object is a full symmetric similarity matrix, not a half-empty one", {
  gli <- testset.gl[1:8, 1:60]
  r <- gl.propShared(gli)

  expect_true(is.matrix(r))
  expect_equal(dim(r), c(8L, 8L))
  # as.matrix(as.dist(res)) mirrors the C++ lower triangle into the upper
  expect_true(isSymmetric(r))
  expect_false(any(r[upper.tri(r)] == 0))
  expect_equal(unname(diag(r)), rep(1, 8))
  expect_equal(rownames(r), indNames(gli))
  expect_equal(colnames(r), indNames(gli))
})

test_that("numeric anchors on testset.gl[1:8, 1:60]", {
  gli <- testset.gl[1:8, 1:60]
  r <- gl.propShared(gli)
  expect_equal(round(r[lower.tri(r)][1:6], 6),
               c(0.981481, 1, 1, 1, 0.990566, 1))
})

test_that("statistic is 1 - mean(|g_i - g_j|)/2 with pairwise-complete NA omission", {
  gli <- testset.gl[1:8, 1:60]
  mat <- as.matrix(gli)
  hand <- matrix(NA_real_, 8, 8)
  for (i in 1:8) {
    for (j in 1:8) {
      d <- abs(mat[i, ] - mat[j, ])
      hand[i, j] <- 1 - mean(d[!is.na(d)]) / 2
    }
  }
  r <- gl.propShared(gli)
  expect_equal(r[lower.tri(r)], hand[lower.tri(hand)])
  # self-comparison is 1 by construction as well as by the diag() override
  expect_equal(unname(diag(hand)), rep(1, 8))
})

test_that("pairs with no overlapping calls return NA, not an error or 1", {
  # [fixed: 1] was NaN from the C++ kernel; gl.dist.ind returns NA. The
  # diagonal is still forced to 1 for the all-NA individual
  set.seed(1)
  m <- matrix(sample(0:2, 6 * 10, replace = TRUE), nrow = 6,
              dimnames = list(paste0("ind", 1:6), paste0("loc", 1:10)))
  m[1, ] <- NA                        # all-NA individual
  m[2, 1:5] <- NA
  m[3, 6:10] <- NA                    # no locus called in both
  gm <- new("genlight", gen = m, ind.names = rownames(m),
            loc.names = colnames(m), ploidy = rep(2, 6))
  gm <- suppressMessages(gl.compliance.check(gm, verbose = 0))

  r <- gl.propShared(gm)
  expect_true(all(is.na(r[1, -1])))
  expect_false(any(is.nan(r)))
  expect_equal(r[1, 1], 1)            # undefined self-similarity reported as 1
  expect_true(is.na(r[2, 3]))
  expect_equal(round(r[4, 5], 6), 0.65)
})

test_that("SilicoDArT input stops with an error", {
  # [fixed: 2] was accepted and returned (1 + simple matching) / 2, a
  # similarity squeezed into [0.5, 1]
  expect_error(gl.propShared(testset.gs[1:6, 1:60], verbose = 0),
               "found SilicoDArT expecting SNP")
})

test_that("non-genlight input fails early with a dartR message", {
  # [fixed: 2] was an S4 dispatch error from indNames() after compiling
  expect_error(gl.propShared(as.matrix(testset.gl[1:3, 1:10]), verbose = 0),
               "inappropriate object passed to function")
})

test_that("verbose controls messaging", {
  # [fixed: 2] the function had no verbose argument and printed nothing
  expect_equal(names(formals(gl.propShared)), c("x", "verbose"))
  out0 <- capture.output(r <- gl.propShared(testset.gl[1:5, 1:50],
                                            verbose = 0))
  expect_length(out0, 0L)
  out1 <- capture.output(r <- gl.propShared(testset.gl[1:5, 1:50],
                                            verbose = 1))
  expect_true(any(grepl("Starting gl.propShared", out1)))
  expect_true(any(grepl("Completed: gl.propShared", out1)))
})

test_that("edge cases: one and two individuals", {
  r1 <- gl.propShared(testset.gl[1, 1:50])
  expect_equal(dim(r1), c(1L, 1L))
  expect_equal(unname(r1[1, 1]), 1)

  r2 <- gl.propShared(testset.gl[1:2, 1:50])
  expect_equal(dim(r2), c(2L, 2L))
  expect_true(isSymmetric(r2))
})

test_that("duplicate individual names produce duplicate dimnames", {
  # [pins defect] name-based lookup silently resolves to the first match
  dup <- testset.gl[1:4, 1:60]
  indNames(dup) <- c("A", "A", "B", "C")
  rd <- gl.propShared(dup)
  expect_equal(colnames(rd), c("A", "A", "B", "C"))
  expect_equal(rd["A", "B"], rd[1, 3])
  expect_false(isTRUE(all.equal(rd[1, 3], rd[2, 3])))
})

test_that("man page is named after the function and joins the distance family", {
  # [fixed: 3] @name was gl.prop.shared and "@family distance" was part of
  # the title
  man <- test_path("..", "..", "man")
  skip_if_not(dir.exists(man), "man/ not available (installed package)")
  expect_false(file.exists(file.path(man, "gl.prop.shared.Rd")))
  expect_true(file.exists(file.path(man, "gl.propShared.Rd")))
  rd <- readLines(file.path(man, "gl.propShared.Rd"), warn = FALSE)
  expect_false(any(grepl("@family", rd, fixed = TRUE)))
  expect_true(any(grepl("\\seealso", rd, fixed = TRUE)))
})

test_that("no C++ is compiled and the result equals 1 - Manhattan distance", {
  # [fixed: 1] the body delegates to gl.dist.ind(method = "manhattan"),
  # exact at every dartR.base release
  expect_false(any(grepl("cppFunction", deparse(body(gl.propShared)))))
  gli <- testset.gl[1:8, 1:60]
  D <- as.matrix(gl.dist.ind(gli, method = "manhattan", type = "matrix",
                             plot.display = FALSE, verbose = 0))
  expect_equal(as.vector(1 - D)[as.vector(lower.tri(D))],
               as.vector(gl.propShared(gli, verbose = 0))[
                 as.vector(lower.tri(D))])
})

test_that("gl.propShared equals 1 - gl.dist.ind(method = 'simple') cell for cell", {
  # Overlap pin. Post-#315 the engine's Simple, Manhattan and Czekanowski
  # are the same statistic; all three complement gl.propShared exactly.
  # Cross-package pin: gl.dist.ind lives in dartR.base. Its "simple" method
  # only equals this statistic once dartR.base PR #315 is merged; skip against
  # an engine that predates it (manhattan/czekanowski match either way).
  .probe <- function() {
    g0 <- testset.gl[1:6, 1:40]
    a <- as.matrix(gl.dist.ind(g0, method = "simple", type = "matrix",
                               plot.display = FALSE, verbose = 0))
    b <- as.matrix(gl.dist.ind(g0, method = "manhattan", type = "matrix",
                               plot.display = FALSE, verbose = 0))
    isTRUE(all.equal(as.vector(a), as.vector(b)))
  }
  skip_if_not(.probe(), "dartR.base engine predates PR #315")
  gli <- testset.gl[1:8, 1:60]
  r <- gl.propShared(gli)
  for (m in c("simple", "manhattan", "czekanowski")) {
    D <- as.matrix(gl.dist.ind(gli, method = m, type = "matrix",
                               plot.display = FALSE, verbose = 0))
    expect_equal(as.vector(1 - D), as.vector(r), info = m)
  }
  # Absolute and Euclidean are different statistics
  Da <- as.matrix(gl.dist.ind(gli, method = "absolute", type = "matrix",
                              plot.display = FALSE, verbose = 0))
  expect_false(isTRUE(all.equal(as.vector(1 - Da), as.vector(r))))
})

test_that("overlap equivalence holds under missing data", {
  # Cross-package pin: gl.dist.ind lives in dartR.base. Its "simple" method
  # only equals this statistic once dartR.base PR #315 is merged; skip against
  # an engine that predates it (manhattan/czekanowski match either way).
  .probe <- function() {
    g0 <- testset.gl[1:6, 1:40]
    a <- as.matrix(gl.dist.ind(g0, method = "simple", type = "matrix",
                               plot.display = FALSE, verbose = 0))
    b <- as.matrix(gl.dist.ind(g0, method = "manhattan", type = "matrix",
                               plot.display = FALSE, verbose = 0))
    isTRUE(all.equal(as.vector(a), as.vector(b)))
  }
  skip_if_not(.probe(), "dartR.base engine predates PR #315")
  set.seed(7)
  m2 <- matrix(sample(0:2, 12 * 40, replace = TRUE), nrow = 12,
               dimnames = list(paste0("i", 1:12), paste0("L", 1:40)))
  m2[sample(length(m2), 120)] <- NA
  gm <- new("genlight", gen = m2, ind.names = rownames(m2),
            loc.names = colnames(m2), ploidy = rep(2, 12))
  gm <- suppressMessages(gl.compliance.check(gm, verbose = 0))

  r <- gl.propShared(gm)
  D <- as.matrix(gl.dist.ind(gm, method = "simple", type = "matrix",
                             plot.display = FALSE, verbose = 0))
  expect_equal(as.vector(1 - D), as.vector(r))
})
