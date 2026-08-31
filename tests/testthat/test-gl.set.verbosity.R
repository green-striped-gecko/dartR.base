# Characterization tests for gl.set.verbosity
# Baseline snapshotted before review (review-gl.set.verbosity).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("valid values set the global option", {
  old <- options(dartR_verbose = NULL); on.exit(options(old))
  o <- capture.output(v <- withVisible(gl.set.verbosity(3)))
  expect_false(v$visible)
  expect_equal(options()$dartR_verbose, 3)
  expect_equal(gl.check.verbosity(), 3)
  # [approved diff F2] baseline: returned NULL although @return
  # promised the value.
  expect_equal(v$value, 3)
  # verbosity 0 sets silently
  o0 <- capture.output(gl.set.verbosity(0))
  expect_length(o0, 0)
  expect_equal(options()$dartR_verbose, 0)
})

test_that("invalid values", {
  old <- options(dartR_verbose = 3); on.exit(options(old))
  # [approved diff F1] baseline: value = 7 warned, silently did NOT set
  # the option, and then printed 'Global verbosity set to: 7' anyway.
  # Now invalid values warn, coerce to 2, and genuinely set 2.
  o <- capture.output(r <- gl.set.verbosity(7))
  expect_equal(options()$dartR_verbose, 2)
  expect_false(any(grepl("set to: 7", o)))
  expect_true(any(grepl("set to: 2", o)))
  expect_equal(r, 2)
  # [approved diff F1] baseline: NULL crashed ('argument is of length
  # zero'); now coerces to 2 like any other invalid value.
  options(dartR_verbose = 3)
  invisible(capture.output(r2 <- gl.set.verbosity(NULL)))
  expect_equal(options()$dartR_verbose, 2)
  expect_equal(r2, 2)
})
