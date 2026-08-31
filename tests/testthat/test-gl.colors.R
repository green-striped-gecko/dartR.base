# Characterization tests for gl.colors
# Baseline snapshotted before review (review-gl.colors).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("fixed palettes return documented colors, silently at verbose 0", {
  o <- capture.output(v <- withVisible(gl.colors(2, verbose = 0)))
  expect_length(o, 0)
  expect_true(v$visible)   # [approved diff F4] the vector is the product
  expect_equal(v$value, c("#3B9AB2", "#78B7C5"))
  expect_equal(gl.colors("2c", verbose = 0), c("deeppink", "chartreuse3"))
  expect_equal(gl.colors("3", verbose = 0),
               c("#3B9AB2", "deeppink", "lemonchiffon"))
  expect_equal(gl.colors(4, verbose = 0),
               c("lemonchiffon", "deeppink", "dodgerblue", "chartreuse3"))
  expect_identical(gl.colors(2, verbose = 0), gl.colors("2", verbose = 0))
})

test_that("palette types return functions with stable output", {
  f.div <- gl.colors("div", verbose = 0)
  f.dis <- gl.colors("dis", verbose = 0)
  f.con <- gl.colors("con", verbose = 0)
  f.vir <- gl.colors("vir", verbose = 0)
  expect_true(all(vapply(list(f.div, f.dis, f.con, f.vir), is.function, TRUE)))
  expect_equal(f.div(4), c("#3B9AB2", "#9EBE91", "#E4B80E", "#F21A00"))
  expect_equal(f.dis(4), c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"))
  expect_equal(f.con(4), c("#0000FF", "#00ABFF", "#FFAC00", "#FF0000"))
  expect_equal(f.vir(4), c("#440154", "#30688D", "#37B778", "#FDE725"))
  expect_length(f.dis(1), 1)
})

test_that("structure palette returns 35 colors", {
  s <- gl.colors("structure", verbose = 0)
  expect_length(s, 35)
  expect_equal(s[1], "#E4E1E3")
})

test_that("invalid types error", {
  expect_error(suppressMessages(capture.output(
    gl.colors("bogus", verbose = 0))))
  expect_error(suppressMessages(capture.output(
    gl.colors(7, verbose = 0))))
  expect_error(suppressMessages(capture.output(
    gl.colors(1, verbose = 0))))
  # BASELINE: documented "pal" type is not accepted
  expect_error(suppressMessages(capture.output(
    gl.colors("pal", verbose = 0))))
})

test_that("error mechanics", {
  # [approved diff F1] baseline: cat(error()) + stop(-1) - the condition
  # message was "-1" and the real message leaked to stdout even inside
  # try(silent = TRUE). Now stop(error()) carries the real message.
  e <- tryCatch(capture.output(gl.colors("bogus", verbose = 0)),
                error = function(e) conditionMessage(e))
  expect_match(e, "No valid color option")
  o <- capture.output(r <- try(gl.colors("bogus", verbose = 0), silent = TRUE))
  expect_length(o, 0)   # nothing to stdout at verbose 0
})

test_that("default-argument evaluation prints banners", {
  # A caller with signature default gl.colors(2) leaks 3 lines when the
  # default is evaluated, regardless of the caller's verbose. Changing
  # the verbose default to 0 was proposed (F2) and REJECTED - the leak
  # is retained deliberately; callers must pass verbose = 0.
  o <- capture.output(x <- gl.colors(2))
  expect_length(o, 3)
})
