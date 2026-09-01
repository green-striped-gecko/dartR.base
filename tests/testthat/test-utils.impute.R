# Characterization tests for the utils.impute helpers
# Baseline snapshotted before review (infrastructure wave, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.

test_that("matrix2gen builds SNPbin lists (serial)", {
  m <- matrix(sample(0:2, 20, replace = TRUE), nrow = 4)
  g <- matrix2gen(m, parallel = FALSE)
  expect_equal(length(g), 4)
  expect_s4_class(g[[1]], "SNPbin")
})

test_that("matrix2gen parallel branch works", {
  m <- matrix(sample(0:2, 20, replace = TRUE), nrow = 4)
  # [approved diff I7] baseline: the parallel branch assigned to i@gen
  # where no object i exists - "object 'i' not found". gl.impute exposes
  # parallel to end users.
  expect_error(matrix2gen(m, parallel = TRUE),
    NA)  # [approved diff I7]
})

test_that("sampling helpers respect their probability contracts", {
  set.seed(1)
  expect_true(is.na(s_alleles(NA)))
  expect_true(is.na(sample_genotype(q_freq = NA)))
  expect_equal(s_alleles(0), 0)   # q = 0: always homozygous reference
  expect_equal(s_alleles(1), 2)   # q = 1: always homozygous alternate
  expect_equal(sample_genotype(q_freq = 0), 0)
  expect_equal(sample_genotype(q_freq = 1), 2)
})
