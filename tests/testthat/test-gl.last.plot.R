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

# ---- base R graphics: recorded from the device ------------------------------

# run code with a graphics device that records (a pdf device keeps no
# display list unless it is enabled); nothing is written to disk
with_recording_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grDevices::dev.control("enable")
  force(code)
}

test_that("a file device without display list records nothing", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(1:3)
  expect_null(utils.plot.record())
})

test_that("a base R plot is kept as a recorded plot; return unchanged", {
  gl.last.plot(clear = TRUE, verbose = 0)
  # at verbose = 0 these functions draw nothing (plot.display is switched off)
  invisible(utils::capture.output(
    res <- with_recording_device(gl.select.shapes(select = 1:5, verbose = 1))))
  expect_equal(res, 1:5)
  p <- gl.last.plot(verbose = 0)
  expect_s3_class(p, "recordedplot")
  # it can be redrawn
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(print(p))
})

test_that("a tree is kept; the tree is still returned", {
  gl.last.plot(clear = TRUE, verbose = 0)
  invisible(utils::capture.output(
    tree <- with_recording_device(gl.tree.nj(testset.gl[1:10, 1:50], verbose = 1))))
  expect_s3_class(tree, "phylo")
  expect_s3_class(gl.last.plot(verbose = 0), "recordedplot")
})

test_that("a heatmap is kept; nothing kept with plot.out = FALSE", {
  D <- dist(as.matrix(testset.gl[1:10, 1:50]))
  invisible(utils::capture.output(with_recording_device(
    suppressWarnings(gl.plot.heatmap(D, verbose = 0)))))
  expect_s3_class(gl.last.plot(verbose = 0), "recordedplot")
  invisible(utils::capture.output(with_recording_device(
    suppressWarnings(gl.plot.heatmap(D, plot.out = FALSE, verbose = 0)))))
  expect_error(gl.last.plot(verbose = 0), "did not draw one")
})

test_that("no plot drawn: an older plot is not returned", {
  x <- testset.gl[1:30, 1:50]
  invisible(gl.report.callrate(x, plot.display = FALSE, verbose = 0))
  expect_true(inherits(gl.last.plot(verbose = 0), "patchwork"))
  invisible(gl.select.shapes(select = 1:5, plot.display = FALSE, verbose = 0))
  expect_error(gl.last.plot(verbose = 0), "gl.select.shapes(), did not draw one",
               fixed = TRUE)
  # a device that cannot record: also no older plot
  invisible(gl.report.callrate(x, plot.display = FALSE, verbose = 0))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  invisible(utils::capture.output(gl.select.shapes(select = 1:5, verbose = 1)))
  expect_error(gl.last.plot(verbose = 0), "did not draw one")
})

test_that("the kept plot is slimmed: same plot, without the function's objects", {
  gl.last.plot(clear = TRUE, verbose = 0)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  # a plotting function whose working space holds a large object (8 MB) that
  # the plot does not use, a label function the plot does use, and a plot
  # that is printed (ggplot2 then caches the mapping in the layers)
  f <- function() {
    big <- stats::rnorm(1e6)
    lab_fun <- function(v) paste0(v, " x")
    d <- data.frame(a = 1:3, b = c(2, 5, 3), g = c("A", "B", "A"))
    p1 <- ggplot2::ggplot(d, ggplot2::aes(a, b, colour = g)) +
      ggplot2::geom_point() +
      ggplot2::scale_x_continuous(labels = lab_fun)
    p2 <- ggplot2::ggplot(d, ggplot2::aes(a)) + ggplot2::geom_histogram(bins = 3)
    p <- patchwork::wrap_plots(p1, p2, ncol = 1)
    print(p)
    utils.plot.keep(p, "plot")
    p
  }
  set.seed(1)
  returned <- f()
  kept <- gl.last.plot(verbose = 0)
  size <- function(o) length(serialize(o, NULL)) / 1e6
  # the returned plot is unchanged and still holds the working space with
  # the large object; the kept plot does not, so it is smaller by at least
  # the size of that object. Compared as a difference, since the size of
  # ggplot2's own objects depends on the ggplot2 and S7 versions.
  big_size <- size(numeric(1e6))
  expect_gt(size(returned), big_size)
  expect_gt(size(returned) - size(kept), big_size)
  # it draws the same: data, and axis labels made by the label function
  for (i in 1:2) {
    expect_equal(ggplot2::layer_data(kept[[i]]), ggplot2::layer_data(returned[[i]]))
  }
  x_labels <- function(g) {
    ggplot2::ggplot_build(g)$layout$panel_params[[1]]$x$get_labels()
  }
  expect_identical(x_labels(kept[[1]]), x_labels(returned[[1]]))
  # and can be changed like any ggplot
  changed <- kept[[1]] + ggplot2::labs(title = "new")
  expect_identical(ggplot2::ggplot_build(changed)$plot$labels$title, "new")
  expect_no_error(print(kept))
})
