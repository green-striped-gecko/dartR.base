# Characterization tests for gl.filter.locmetric
# Baseline snapshotted before review (review-gl.filter.locmetric).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

make_test <- function() {
  test <- testset.gl
  test@other$loc.metrics$test <- 1:nLoc(test)
  test
}

test_that("gl.filter.locmetric keep='within' retains [lower, upper] inclusive", {
  test <- make_test()
  res <- gl.filter.locmetric(test, metric = "test", upper = 255, lower = 200,
                             keep = "within", verbose = 0)
  expect_equal(nLoc(res), 56)
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
  expect_true(all(res@other$loc.metrics$test >= 200 &
                    res@other$loc.metrics$test <= 255))
  expect_true(200 %in% res@other$loc.metrics$test)  # boundary kept
})

test_that("gl.filter.locmetric appends history and leaves input untouched", {
  test <- make_test()
  h <- length(test@other$history)
  res <- gl.filter.locmetric(test, metric = "test", upper = 255, lower = 200,
                             verbose = 0)
  expect_equal(length(res@other$history), h + 1)
  test2 <- make_test()
  invisible(gl.filter.locmetric(test2, metric = "test", upper = 255,
                                lower = 200, verbose = 0))
  expect_identical(test2, make_test())
})

test_that("gl.filter.locmetric errors when the metric is absent", {
  expect_error(
    gl.filter.locmetric(testset.gl, metric = "nope", upper = 1, lower = 0,
                        verbose = 0),
    "not found"
  )
})

test_that("gl.filter.locmetric with NA metric values (within)", {
  # NA-metric loci are removed (which() drops NA comparisons); sync holds.
  tna <- make_test()
  tna@other$loc.metrics$test[201:210] <- NA
  res <- gl.filter.locmetric(tna, metric = "test", upper = 255, lower = 200,
                             keep = "within", verbose = 0)
  expect_equal(nLoc(res), 46)
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
  expect_equal(sum(is.na(res@other$loc.metrics$test)), 0)
})

test_that("gl.filter.locmetric keep='outside'", {
  # [approved diff F1] Pre-fix, 'outside' ALWAYS crashed (impossible AND
  # condition). Post-fix it retains the exact complement of 'within':
  # metric < lower or > upper.
  test <- make_test()
  res <- gl.filter.locmetric(test, metric = "test", upper = 255, lower = 200,
                             keep = "outside", verbose = 0)
  expect_equal(nLoc(res), sum(1:255 < 200))
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
  expect_true(all(res@other$loc.metrics$test < 200 |
                    res@other$loc.metrics$test > 255))
  # 'within' and 'outside' partition the loci exactly
  res_w <- gl.filter.locmetric(test, metric = "test", upper = 255,
                               lower = 200, keep = "within", verbose = 0)
  expect_equal(nLoc(res) + nLoc(res_w), nLoc(test))
})

test_that("gl.filter.locmetric with an invalid keep value", {
  # [approved diff F2] invalid keep now coerces to 'within' with a gated
  # warning instead of crashing with "object 'x2' not found".
  test <- make_test()
  out <- capture.output(
    res <- gl.filter.locmetric(test, metric = "test", upper = 255,
                               lower = 200, keep = "bogus", verbose = 0)
  )
  expect_length(out, 0)
  expect_equal(nLoc(res), 56)
})

test_that("gl.filter.locmetric with a non-numeric metric", {
  # [approved diff F2] non-numeric metrics now stop with a clear fatal
  # error instead of factor warnings + a zero-loci crash.
  test <- make_test()
  expect_error(
    gl.filter.locmetric(test, metric = "AlleleID", upper = 255,
                        lower = 200, verbose = 0),
    "not numeric"
  )
})

test_that("gl.filter.locmetric pop side effect on pop-less objects", {
  # [approved diff F3] pop-less objects keep their pop state; the 'pop1'
  # stamping preamble is removed.
  test <- make_test()
  test@pop <- NULL
  res <- gl.filter.locmetric(test, metric = "test", upper = 255, lower = 200,
                             verbose = 0)
  expect_true(is.null(res@pop) || length(res@pop) == 0)
})

test_that("gl.filter.locmetric returns invisibly", {
  # [approved diff F5] now returns invisibly, matching the other filters.
  test <- make_test()
  v <- withVisible(gl.filter.locmetric(test, metric = "test", upper = 255,
                                       lower = 200, verbose = 0))
  expect_false(v$visible)
})
