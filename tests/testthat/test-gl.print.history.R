# Characterization tests for gl.print.history
# Baseline snapshotted before review (review-gl.print.history).
# Assertions marked [approved diff] were flipped in Phase C.

gl3 <- gl.filter.callrate(testset.gl, method = "loc", threshold = 0.9,
                          verbose = 0)
gl3 <- gl.filter.callrate(gl3, method = "ind", threshold = 0.95,
                          verbose = 0)

test_that("full history prints one numbered entry per call", {
  o <- capture.output(v <- withVisible(gl.print.history(gl3)))
  expect_false(v$visible)
  # [approved diff changes 1, 3] baseline: knitr markdown table, with the
  # wrapped first call spilling onto a bare line "probar = TRUE) |".
  expect_true(any(grepl("^1 gl.read.dart", o)))
  expect_true(any(grepl("^3 gl.filter.callrate", o)))
  expect_true(any(grepl("^    probar = TRUE\\)$", o)))
  expect_false(any(grepl("|", o, fixed = TRUE)))
})

test_that("the history table is returned at every verbosity", {
  # [approved diff change 5] baseline: returned NULL; verbose = 0 did
  # nothing at all.
  o <- capture.output(v <- gl.print.history(gl3, verbose = 0))
  expect_length(o, 0)
  expect_s3_class(v, "data.frame")
  expect_equal(v$nr, 1:3)
  expect_match(v$history[3], "^gl.filter.callrate\\(x = gl3")
})

test_that("selected entries keep their numbers; out of range errors", {
  # [approved diff change 2] baseline: history = 3 was labelled 1, and
  # c(1, 7) printed a "2 | NULL" row.
  o <- capture.output(gl.print.history(gl3, history = 3))
  expect_true(any(grepl("^3 gl.filter.callrate", o)))
  expect_error(gl.print.history(gl3, history = c(1, 7), verbose = 0),
               "entry numbers from 1 to 3")
})

test_that("a history list is printed with or without x", {
  o <- capture.output(gl.print.history(history = gl3@other$history))
  expect_true(any(grepl("^3 gl.filter.callrate", o)))
  # [approved diff change 4] baseline: with x given, a list stopped with
  # "invalid subscript type 'list'".
  v <- gl.print.history(gl3, history = gl3@other$history, verbose = 0)
  expect_equal(nrow(v), 3)
})

test_that("invalid inputs stop with a clear message", {
  # [approved diff change 4] baseline: "object 'hist2' not found"
  expect_error(gl.print.history(verbose = 0), "provide a genlight object")
  expect_error(gl.print.history(data.frame(a = 1), verbose = 0),
               "provide a genlight object")
})

test_that("an empty history warns and returns an empty table", {
  o <- capture.output(v <- gl.print.history(history = list(), verbose = 2))
  expect_true(any(grepl("no history entries found", o)))
  expect_equal(nrow(v), 0)
})
