# Characterization tests for utils.stats (std.error)
# Baseline snapshotted before review (infrastructure wave, dev at ddaed27).

test_that("std.error matches the definition", {
  x <- c(1, 2, 3, NA)
  expect_equal(std.error(x), sqrt(var(x, na.rm = TRUE) / 3))
  expect_equal(std.error(c(5, 5, 5)), 0)
})
