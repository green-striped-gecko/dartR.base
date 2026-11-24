# gl.last.plot() and utils.plot.keep(): the plots of the last plotting
# function can be recovered without saving them, for both storage backends,
# without changing what the plotting functions return.

test_that("no plot yet, or cleared, gives an informative error", {
  gl.last.plot(clear = TRUE, verbose = 0)
  expect_error(gl.last.plot(verbose = 0), "no dartR function has made a plot")
})

test_that("objects kept in one call form one set; the next call replaces it", {
  gl.last.plot(clear = TRUE, verbose = 0)
  two <- function() {
    utils.plot.keep(ggplot2::ggplot(), "plot")
    utils.plot.keep(data.frame(a = 1), "table")
    invisible(NULL)
  }
  one <- function() {
    utils.plot.keep(ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)), "plot")
    invisible(NULL)
  }
  two()
  res <- gl.last.plot(verbose = 0)
  expect_type(res, "list")
  expect_named(res, c("plot", "table"))
  expect_s3_class(gl.last.plot("table", verbose = 0), "data.frame")
  expect_true(inherits(gl.last.plot(1, verbose = 0), "ggplot"))
  expect_error(gl.last.plot("nope", verbose = 0), "which must be one of")
  expect_error(gl.last.plot(3, verbose = 0), "which must be one of")

  one()
  res <- gl.last.plot(verbose = 0)
  expect_true(inherits(res, "ggplot"))
  # an aesthetic is a quosure, i.e. a one-sided formula
  expect_equal(deparse(unclass(res$mapping$x)[[2]]), "mpg")

  # a second call of the same function starts a new set
  two(); two()
  expect_length(gl.last.plot(verbose = 0), 2)

  expect_null(gl.last.plot(clear = TRUE, verbose = 0))
  expect_error(gl.last.plot(verbose = 0), "no dartR function")
})

test_that("invalid arguments are refused", {
  expect_error(gl.last.plot(clear = NA, verbose = 0), "clear must be TRUE or FALSE")
  expect_error(utils.plot.keep(1, name = ""), "name must be")
})

test_that("a report keeps its plot without plot.display or plot.file; return unchanged", {
  gl.last.plot(clear = TRUE, verbose = 0)
  x <- testset.gl[1:30, 1:50]
  out <- gl.report.callrate(x, plot.display = FALSE, verbose = 0)
  expect_identical(out, x)
  p <- gl.last.plot(verbose = 0)
  expect_true(inherits(p, "patchwork"))
  expect_equal(length(p), 2)
})

test_that("a filter returns the filtered object and keeps its plot", {
  gl.last.plot(clear = TRUE, verbose = 0)
  x <- testset.gl[1:30, 1:50]
  out <- gl.filter.callrate(x, threshold = 0.9, plot.display = TRUE,
                            verbose = 0)
  expect_s4_class(out, "genlight")
  expect_true(nLoc(out) <= nLoc(x))
  expect_true(inherits(gl.last.plot(verbose = 0), "patchwork"))
})

test_that("the kept plot is the plot that plot.file saves", {
  gl.last.plot(clear = TRUE, verbose = 0)
  x <- testset.gl[1:30, 1:50]
  d <- tempfile("lastplot")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  invisible(gl.report.callrate(x, plot.display = FALSE, plot.dir = d,
                               plot.file = "cr", verbose = 0))
  saved <- readRDS(file.path(d, "cr.RDS"))
  kept <- gl.last.plot(verbose = 0)
  expect_equal(ggplot2::layer_data(kept[[1]]), ggplot2::layer_data(saved[[1]]))
  expect_equal(ggplot2::layer_data(kept[[2]]), ggplot2::layer_data(saved[[2]]))
})

test_that("both storage backends keep the same plot", {
  x <- testset.gl[1:30, 1:50]
  invisible(gl.report.callrate(x, plot.display = FALSE, verbose = 0))
  p_gen <- gl.last.plot(verbose = 0)
  xf <- gl.gen2fbm(x, verbose = 0)
  invisible(gl.report.callrate(xf, plot.display = FALSE, verbose = 0))
  p_fbm <- gl.last.plot(verbose = 0)
  expect_equal(ggplot2::layer_data(p_fbm[[1]]), ggplot2::layer_data(p_gen[[1]]))
  expect_equal(ggplot2::layer_data(p_fbm[[2]]), ggplot2::layer_data(p_gen[[2]]))
})

test_that("a function that makes a plot and a table returns both", {
  gl.last.plot(clear = TRUE, verbose = 0)
  x <- testset.gl[1:30, 1:50]
  # the function prints its plot; ggplot2 warns about points outside the
  # axis range, which is not what is tested here
  invisible(utils::capture.output(suppressWarnings(
    gl.diagnostics.hwe(x, verbose = 0))))
  res <- gl.last.plot(verbose = 0)
  expect_named(res, c("plot", "table"))
  expect_true(inherits(res$plot, "ggplot") || inherits(res$plot, "patchwork"))
  expect_s3_class(res$table, "data.frame")
})
