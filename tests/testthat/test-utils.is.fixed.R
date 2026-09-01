# Characterization tests for utils.is.fixed
# Baseline snapshotted before review (kernel-wave review, dev at ddaed27).

test_that("truth table and numeric return contract", {
  expect_identical(utils.is.fixed(100, 0), 1)
  expect_identical(utils.is.fixed(0, 100), 1)
  expect_identical(utils.is.fixed(80, 0), 0)
  expect_identical(utils.is.fixed(50, 50), 0)
  expect_true(is.na(utils.is.fixed(NA, 0)))
  expect_true(is.na(utils.is.fixed(100, NA)))
  # tolerance boundaries (docs: 95,5 at tloc = 0.05 is fixed)
  expect_identical(utils.is.fixed(95, 5, tloc = 0.05), 1)
  expect_identical(utils.is.fixed(96, 4, tloc = 0.05), 1)
  expect_identical(utils.is.fixed(94, 6, tloc = 0.05), 0)
  # the return is numeric 1/0 (callers test > 0), not logical
  expect_type(utils.is.fixed(100, 0), "double")
})
