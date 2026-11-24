# gl.plot.help(): the code it prints refers to the plot's own labels, groups
# and columns, and it parses and runs.

# run the suggested code (except saving) with the plot bound to its name
run_code <- function(code, ...) {
  env <- list2env(list(...), parent = asNamespace("dartR.base"))
  code <- code[!grepl("^#|^\\s*$|ggsave\\(|gl.last.plot\\(", code)]
  for (line in code) {
    res <- eval(parse(text = line), envir = env)
    if (inherits(res, "ggplot")) {
      # building the plot shows that the added layers and scales work
      suppressMessages(suppressWarnings(ggplot2::ggplot_build(res)))
    }
  }
  invisible(TRUE)
}

test_that("code for a ggplot names its labels and groups, parses and runs", {
  d <- data.frame(PC1 = c(1, 2, 3, 4), PC2 = c(2, 1, 4, 3),
                  pop = c("A", "A", "B", "B"), ind = c("i1", "i2", "i3", "i4"))
  p <- ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, colour = pop, shape = pop)) +
    ggplot2::geom_point() + ggplot2::labs(title = "PCA")
  code <- gl.plot.help(p, verbose = 0)
  expect_type(code, "character")
  expect_true(any(grepl('labs(x = "PC1", y = "PC2", title = "PCA")', code, fixed = TRUE)))
  expect_true(any(grepl('scale_colour_manual(values = c("A" = ', code, fixed = TRUE)))
  expect_true(any(grepl("scale_shape_manual", code, fixed = TRUE)))
  expect_true(any(grepl('subset(p$data, ind == "i1")', code, fixed = TRUE)))
  expect_true(any(grepl("ggsave(", code, fixed = TRUE)))
  expect_silent(parse(text = code))
  expect_true(run_code(code, p = p))
})

test_that("code for a patchwork addresses each panel and runs", {
  a <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg)) +
    ggplot2::geom_histogram(bins = 5, fill = "blue")
  b <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  fig <- a / b
  code <- gl.plot.help(fig, verbose = 0)
  expect_true(any(grepl("fig & theme(", code, fixed = TRUE)))
  expect_true(any(grepl("fig[[1]] <- fig[[1]] + labs(", code, fixed = TRUE)))
  expect_true(any(grepl("fig[[2]] <- fig[[2]] + labs(", code, fixed = TRUE)))
  expect_true(any(grepl('fig[[1]]$layers[[1]]$aes_params$fill <- "blue"', code, fixed = TRUE)))
  expect_silent(parse(text = code))
  expect_true(run_code(code, fig = fig))
})

test_that("more than 10 groups give valid code", {
  d <- data.frame(x = 1:12, y = 1:12, g = LETTERS[1:12])
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y, colour = g)) + ggplot2::geom_point()
  code <- gl.plot.help(p, verbose = 0)
  expect_true(any(grepl("first 10 of 12 groups", code, fixed = TRUE)))
  expect_silent(parse(text = code))
})

test_that("a list from gl.last.plot() skips tables; default is the last plot", {
  gl.last.plot(clear = TRUE, verbose = 0)
  f <- function() {
    utils.plot.keep(ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
                      ggplot2::geom_point(), "plot")
    utils.plot.keep(data.frame(a = 1), "table")
  }
  f()
  code <- gl.plot.help(verbose = 0)
  expect_equal(code[2], "p <- gl.last.plot()")
  expect_true(any(grepl("p$plot + labs(", code, fixed = TRUE)))
  expect_false(any(grepl("p$table", code, fixed = TRUE)))
  expect_silent(parse(text = code))
})

test_that("a dartR plot gives runnable code", {
  x <- testset.gl[1:30, 1:50]
  invisible(gl.report.callrate(x, plot.display = FALSE, verbose = 0))
  p <- gl.last.plot(verbose = 0)
  code <- gl.plot.help(p, verbose = 0)
  expect_true(run_code(code, p = p))
})

test_that("printing follows verbose; invalid input is refused", {
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  expect_silent(gl.plot.help(p, verbose = 0))
  expect_output(gl.plot.help(p, verbose = 1), "labs\\(x = \"wt\"")
  expect_error(gl.plot.help(1:3, verbose = 0), "p must be a ggplot")
  expect_error(gl.plot.help(list(a = 1), verbose = 0), "contains no ggplot")
})
