# .as_dartR() (utils.dartR.class.def.r) must return a valid dartR object for
# every input; `class(x) <- "dartR"` gave objects without the fbm slot that
# failed validObject().
plain_genlight <- function() {
  g <- new("genlight", as.matrix(platypus.gl), ploidy = 2)
  pop(g) <- pop(platypus.gl)
  g@other <- platypus.gl@other
  g
}
relabelled <- function() {
  x <- plain_genlight()
  class(x) <- "dartR"
  x
}

test_that("a plain genlight is coerced to a valid dartR with every slot kept", {
  g <- plain_genlight()
  x <- .as_dartR(g)
  expect_s4_class(x, "dartR")
  expect_true(methods::.hasSlot(x, "fbm"))
  expect_true(validObject(x))
  for (s in slotNames("genlight")) {
    expect_identical(slot(x, s), slot(g, s), label = s)
  }
})

test_that("a relabelled object without the fbm slot is rebuilt", {
  b <- relabelled()
  expect_false(methods::.hasSlot(b, "fbm"))
  expect_error(validObject(b), "fbm")
  x <- .as_dartR(b)
  expect_true(methods::.hasSlot(x, "fbm"))
  expect_true(validObject(x))
  expect_identical(as.matrix(x), as.matrix(plain_genlight()))
  expect_identical(x@other, platypus.gl@other)
})

test_that("packaged datasets are rebuilt; a valid dartR object is returned unchanged", {
  # dartR.data objects were saved before the fbm slot existed
  expect_false(methods::.hasSlot(platypus.gl, "fbm"))
  x <- .as_dartR(platypus.gl)
  expect_true(validObject(x))
  for (s in slotNames("genlight")) {
    expect_identical(slot(x, s), slot(platypus.gl, s), label = s)
  }
  expect_identical(.as_dartR(x), x)
  skip_if_not_installed("bigstatsr")
  capture.output(f <- gl.gen2fbm(testset.gl, verbose = 0))
  expect_identical(.as_dartR(f), f)
})

test_that("gl.compliance.check returns a valid object for both inputs", {
  capture.output(a <- gl.compliance.check(plain_genlight(), verbose = 0))
  expect_true(validObject(a))
  capture.output(b <- gl.compliance.check(relabelled(), verbose = 0))
  expect_true(validObject(b))
})

test_that("gl.load repairs an object saved without the fbm slot", {
  f <- tempfile(fileext = ".rds")
  saveRDS(relabelled(), f)
  capture.output(x <- gl.load(f, verbose = 0))
  expect_true(validObject(x))
  expect_identical(as.matrix(x), as.matrix(plain_genlight()))
})

test_that("functions that converted with class<- now return valid objects", {
  g <- plain_genlight()
  # no dartR metadata, so that adegenet's `[` keeps the object consistent
  g0 <- new("genlight", as.matrix(platypus.gl)[, 1:20], ploidy = 2)
  pop(g0) <- pop(platypus.gl)
  capture.output({
    k <- gl.keep.pop(g, pop.list = popNames(g)[1], verbose = 0)
    d <- gl.drop.ind(g, ind.list = indNames(g)[1:3], verbose = 0)
    s <- gl.sort(g, verbose = 0)
    j <- gl.join(g0[, 1:10], g0[, 11:20], verbose = 0)
  })
  for (o in list(k, d, s, j)) expect_true(validObject(o))
})
