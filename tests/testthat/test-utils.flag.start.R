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

test_that("the function name is printable for pkg::fun and do.call(fun)", {
  # match.call()[[1]] is a name, a call (pkg::fun) or the function itself
  expect_identical(utils.flag.start(as.name("gl.x"), verbose = 0), "gl.x")
  expect_identical(utils.flag.start(quote(dartR.base::gl.sort), verbose = 0),
                   "gl.sort")
  expect_identical(utils.flag.start(gl.filter.callrate, verbose = 0),
                   "gl.filter.callrate")
  o <- capture.output(utils.flag.start(gl.filter.callrate, verbose = 1))
  expect_true(grepl("Starting gl.filter.callrate", o[1]))
  # a whole run through do.call prints (it stopped at cat() before)
  o <- capture.output(do.call(gl.filter.callrate,
                              list(x = testset.gl, threshold = 0.9,
                                   verbose = 1)))
  expect_true(any(grepl("Completed: gl.filter.callrate", o)))
  o <- capture.output(dartR.base::gl.filter.callrate(testset.gl,
                                                     threshold = 0.9,
                                                     verbose = 1))
  expect_identical(o[length(o)], "Completed: gl.filter.callrate ")
})

test_that("utils.history.call keeps names and small values, not objects", {
  gl <- testset.gl
  a <- gl.filter.callrate(gl, threshold = 0.9, verbose = 0)
  ha <- a@other$history[[length(a@other$history)]]
  expect_identical(ha, quote(gl.filter.callrate(x = gl, threshold = 0.9,
                                                verbose = 0)))
  b <- do.call(gl.filter.callrate,
               list(x = gl, threshold = 0.9, verbose = 0))
  hb <- b@other$history[[length(b@other$history)]]
  expect_identical(hb[[1]], as.name("gl.filter.callrate"))
  expect_true(is.name(hb$x))
  expect_identical(as.character(hb$x), "<dartR: 250 ind x 255 loc>")
  expect_identical(hb$threshold, 0.9)
  expect_identical(as.matrix(a), as.matrix(b))
  # the history no longer grows with the object under chained do.call
  x <- gl
  for (i in 1:3) {
    x <- do.call(gl.filter.callrate, list(x = x, threshold = 0.5, verbose = 0))
  }
  expect_lt(as.numeric(object.size(x@other$history)), 20000)
  # other values: labels; NULL and expressions stay
  h <- utils.history.call(as.call(list(as.name("f"), d = data.frame(a = 1:3),
                                       v = 1:100, n = NULL, e = quote(y + 1))))
  expect_identical(as.character(h$d), "<data.frame: 3 x 1>")
  expect_identical(as.character(h$v), "<integer of length 100>")
  expect_identical(h$e, quote(y + 1))
  expect_true("n" %in% names(h))
})
