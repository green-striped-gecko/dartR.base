# Characterization tests for utils.flag.start
# Baseline snapshotted before review (infrastructure wave, dev at ddaed27).

test_that("start flag prints per verbosity contract", {
  expect_equal(length(capture.output(utils.flag.start("gl.x", verbose = 0))), 0)
  o1 <- capture.output(utils.flag.start("gl.x", verbose = 1))
  expect_true(grepl("Starting gl.x", o1[1]))
  o5 <- capture.output(utils.flag.start("gl.x", build = "B", verbose = 5))
  expect_true(any(grepl("Build", o5)))
  expect_error(utils.flag.start(NULL, verbose = 0), "calling function")
  r <- withVisible(utils.flag.start("gl.x", verbose = 0))
  expect_false(r$visible)
  expect_equal(r$value, "gl.x")
})
