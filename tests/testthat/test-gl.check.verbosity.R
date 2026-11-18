# Characterization tests for gl.check.verbosity
# Baseline snapshotted before review (review-gl.check.verbosity).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("NULL falls back to the global option, then to 2", {
  old <- options(dartR_verbose = NULL); on.exit(options(old))
  expect_equal(gl.check.verbosity(), 2)
  expect_equal(gl.check.verbosity(NULL), 2)
  options(dartR_verbose = 3)
  expect_equal(gl.check.verbosity(), 3)
  options(dartR_verbose = 0)
  expect_equal(gl.check.verbosity(), 0)
})

test_that("an explicit value overrides the global option", {
  old <- options(dartR_verbose = 5); on.exit(options(old))
  expect_equal(gl.check.verbosity(1), 1)
  expect_equal(gl.check.verbosity(0), 0)
})

test_that("values 0 to 5 pass through silently", {
  for (v in 0:5) {
    o <- capture.output(r <- gl.check.verbosity(v))
    expect_equal(r, v)
    expect_length(o, 0)
  }
  # non-integers inside the range are accepted unchanged
  expect_equal(gl.check.verbosity(2.5), 2.5)
})

test_that("out-of-range and non-numeric values warn and return 2", {
  for (v in list(6, -1, TRUE, "3", NA)) {
    o <- capture.output(r <- gl.check.verbosity(v))
    expect_equal(r, 2)
    # [approved diff change 3] baseline: "must be an integer in the
    # range \n<20 spaces>0 to 5" split over two lines.
    expect_length(o, 1)
    expect_true(grepl("must be a single number from 0 to 5", o))
  }
  o <- capture.output(r <- gl.check.verbosity(6))
  expect_true(grepl("(received 6); set to 2", o, fixed = TRUE))
})

test_that("NA_real_ and wrong-length values warn and return 2", {
  # [approved diff change 1] baseline: these stopped with R-internal
  # errors ("missing value where TRUE/FALSE needed", "the condition has
  # length > 1", "argument is of length zero").
  for (v in list(NA_real_, c(1, 3), numeric(0))) {
    o <- capture.output(r <- gl.check.verbosity(v))
    expect_equal(r, 2)
    expect_true(grepl("verbose must be a single number", o))
  }
})

test_that("an invalid global option warns and returns 2", {
  # [approved diff change 2] baseline: "loud" and 9 were returned as is.
  old <- options(dartR_verbose = "loud"); on.exit(options(old))
  o <- capture.output(r <- gl.check.verbosity())
  expect_equal(r, 2)
  expect_true(grepl("option dartR_verbose must be", o))
  options(dartR_verbose = 9)
  o <- capture.output(r <- gl.check.verbosity())
  expect_equal(r, 2)
  # an explicit valid argument still wins over an invalid option
  o <- capture.output(r <- gl.check.verbosity(1))
  expect_equal(r, 1)
  expect_length(o, 0)
})

test_that("a caller given verbose = NA_real_ runs instead of stopping", {
  # [approved diff change 1] baseline: stopped inside gl.report.callrate
  invisible(capture.output(
    res <- gl.report.callrate(testset.gl, plot.display = FALSE,
                              verbose = NA_real_)
  ))
  expect_s4_class(res, "genlight")
})
