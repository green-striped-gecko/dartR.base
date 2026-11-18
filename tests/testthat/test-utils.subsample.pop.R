# Characterization tests for utils.subsample.pop (internal; used by
# gl.report.heterozygosity and gl.report.polyploid_heterozygosity via
# subsample.pop = TRUE), captured before the function-review changes and
# updated for the approved changes (report:
# function-review/reports/dartR.base/utils.subsample.pop.md).

test_that("default output on testset.gl", {
  set.seed(1)
  a <- utils.subsample.pop(testset.gl, n.limit = 10)
  expect_equal(dim(a), c(100L, 4L))
  expect_named(a, c("res.mean", "res_SE", "pop", "subsample"))
  expect_equal(unique(a$subsample), c(10, 5, 4, 3, 2))
  # only populations with at least n.limit individuals
  expect_setequal(unique(a$pop),
                  names(which(table(pop(testset.gl)) >= 10)))
  expect_equal(a$res.mean[1], 0.01356009, tolerance = 1e-6)
})

test_that("a subsample of the whole population reproduces its Ho", {
  p <- names(which(table(pop(testset.gl)) >= 10))[1]
  y <- testset.gl[pop(testset.gl) == p, ]
  m <- as.matrix(y)
  set.seed(2)
  r <- utils.subsample.pop(y, n.limit = 10, subsamples = nrow(m))
  expect_equal(r$res.mean, mean(colMeans(m == 1, na.rm = TRUE), na.rm = TRUE))
  expect_equal(r$res_SE, 0)
})

test_that("FBM-backed input gives the dense result", {
  skip_if_not_installed("bigstatsr")
  set.seed(1)
  a <- utils.subsample.pop(testset.gl, n.limit = 10)
  set.seed(1)
  b <- utils.subsample.pop(gl.gen2fbm(testset.gl, verbose = 0), n.limit = 10)
  expect_identical(a, b)
})

test_that("edge cases", {
  # change 1: n.limit below the largest subsample size previously errored;
  # each population now uses only the sizes it can supply
  set.seed(1)
  a <- utils.subsample.pop(testset.gl, n.limit = 5)
  tab <- table(pop(testset.gl))
  expect_setequal(unique(a$pop), names(which(tab >= 5)))
  small <- names(which(tab >= 5 & tab < 10))[1]
  expect_equal(a$subsample[a$pop == small], c(5, 4, 3, 2))
  expect_true(all(a$subsample <= tab[a$pop]))
  # the same seed gives the n.limit = 10 result when all sizes fit
  set.seed(1)
  b10 <- utils.subsample.pop(testset.gl, n.limit = 10)
  expect_equal(b10$res.mean[1], 0.01356009, tolerance = 1e-6)
  # via the report
  pdf(NULL)
  on.exit(dev.off())
  r <- gl.report.heterozygosity(testset.gl, subsample.pop = TRUE,
                                n.limit = 5, verbose = 0)
  expect_named(r, c("subsample", "results"))
  # change 2: a single locus previously errored
  set.seed(1)
  expect_equal(nrow(utils.subsample.pop(testset.gl[, 1], n.limit = 10)),
               100L)
  # no population reaches n.limit: empty table
  expect_equal(nrow(utils.subsample.pop(testset.gl, n.limit = 1000)), 0L)
})
