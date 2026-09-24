# Characterization tests for gl.fdsim
# Baseline snapshotted before review (review-gl.fdsim).
# Assertions marked [approved diff] were flipped in Phase C.

pp <- c("EmsubRopeMata", "EmmacBurnBara")

test_that("allopatric output: named list of four scalars", {
  set.seed(1)
  o <- gl.fdsim(testset.gl, poppair = pp, reps = 200, verbose = 0)
  expect_type(o, "list")
  expect_named(o, c("observed", "mnexpected", "sdexpected", "prob"))
  expect_equal(o$observed, 12)
  # [approved diff changes 1, 5] baseline: mnexpected 2.326186,
  # sdexpected 0.8092613 (spread of the expectation), prob 3.098113e-33
  # (normal curve). Now the simulated count and its tail.
  expect_equal(o$mnexpected, 2.355)
  expect_equal(o$sdexpected, 1.202583, tolerance = 1e-6)
  expect_equal(o$prob, 1 / 201)
})

test_that("p-value is the simulated tail of the false-positive count", {
  # [approved diff change 1] baseline: obs = 4 gave p = 0.0218 from
  # pnorm(4, 2.345, 0.820, lower.tail = FALSE).
  set.seed(5)
  o <- gl.fdsim(testset.gl, poppair = pp, obs = 4, reps = 1000,
                verbose = 0)
  expect_equal(o$observed, 4)
  expect_equal(o$prob, 0.1468531, tolerance = 1e-6)
  # an independent re-implementation of the model (5000 counts) gives
  # P(count >= 4) = 0.152; allow for simulation error
  expect_gt(o$prob, 0.12)
  expect_lt(o$prob, 0.18)
})

test_that("rare false positives: p no longer collapses to ~0", {
  # [approved diff change 1] baseline: obs = 1 gave p = 4.6e-20;
  # the simulated count reaches 1 in about 4% of replicates.
  set.seed(5)
  o <- gl.fdsim(platypus.gl, poppair = c("SEVERN_BELOW", "TENTERFIELD"),
                obs = 1, reps = 5000, verbose = 0)
  expect_gt(o$prob, 0.03)
  expect_lt(o$prob, 0.055)
})

test_that("sympatric pairs use the pooled frequency, symmetric in order", {
  # [approved diff change 2] baseline: rows from the interleaved table,
  # mnexpected 0.001559357 and identical output for both orders.
  s <- replicate(2, {
    set.seed(1)
    gl.fdsim(testset.gl, poppair = pp, reps = 2000, sympatric = TRUE,
             verbose = 0)$mnexpected
  })
  set.seed(2)
  s2 <- gl.fdsim(testset.gl, poppair = rev(pp), reps = 2000,
                 sympatric = TRUE, verbose = 0)$mnexpected
  expect_lt(s[1], 0.05)
  expect_lt(abs(s[1] - s2), 0.02)
})

test_that("SilicoDArT sample size counts individuals", {
  # [approved diff change 3] baseline: runs with 2n; now n
  set.seed(1)
  o <- gl.fdsim(testset.gs, poppair = pp, reps = 200, verbose = 0)
  expect_equal(o$observed, 47)
  expect_equal(o$mnexpected, 6.6)
})

test_that("argument errors", {
  # [approved diff change 4] baseline: length-1 poppair said
  # "Population B mislabelled"; c(A, A) and reps = 1 ran; obs unchecked.
  expect_error(gl.fdsim(testset.gl, poppair = "EmsubRopeMata", verbose = 0),
               "two different population labels")
  expect_error(gl.fdsim(testset.gl, poppair = c(pp[1], pp[1]),
                        verbose = 0),
               "two different population labels")
  expect_error(gl.fdsim(testset.gl, poppair = c("nope", pp[2]),
                        verbose = 0),
               "population nope not found")
  expect_error(gl.fdsim(testset.gl, poppair = pp, reps = 1, verbose = 0),
               "reps must be a whole number")
  expect_error(gl.fdsim(testset.gl, poppair = pp, obs = -1, verbose = 0),
               "obs must be NULL")
})

test_that("verbose 3 labels sample sizes with population names", {
  # [approved diff change 7] baseline: "Sample sizes: 11 5", factor order
  set.seed(1)
  o <- capture.output(gl.fdsim(testset.gl[, 1:100], poppair = pp,
                               reps = 20, verbose = 3))
  expect_true(any(grepl("EmsubRopeMata = 5, EmmacBurnBara = 11", o)))
})

test_that("FBM-backed input gives the same result", {
  f <- gl.gen2fbm(testset.gl, verbose = 0)
  set.seed(1)
  a <- gl.fdsim(f, poppair = pp, reps = 50, verbose = 0)
  set.seed(1)
  b <- gl.fdsim(testset.gl, poppair = pp, reps = 50, verbose = 0)
  expect_equal(a, b)
})
