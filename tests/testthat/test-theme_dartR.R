# Characterization tests for theme_dartR, captured before the
# function-review changes and updated for the approved changes (report:
# function-review/reports/dartR.base/theme_dartR.md).

test_that("theme object and key elements", {
  th <- theme_dartR()
  expect_s3_class(th, "theme")
  expect_true(attr(th, "complete"))
  expect_equal(th$text$size, 11)
  expect_equal(th$line$linewidth, 0.5)
  expect_equal(th$axis.ticks.length, grid::unit(2.75, "pt"))
  expect_identical(th$plot.title$face, "bold")
  expect_identical(th$legend.position, "right")
  expect_equal(theme_dartR(base_size = 20)$text$size, 20)
})

test_that("facet strip size", {
  # change 1: previously fixed at 14 pt whatever base_size was; now 14 pt at
  # the default and scaled otherwise
  strip_pt <- function(bs) {
    p <- ggplot2::ggplot(data.frame(a = 1, g = "x"), ggplot2::aes(a, a)) +
      ggplot2::facet_wrap(~g) + theme_dartR(base_size = bs)
    ggplot2::calc_element("strip.text.x.top",
                          ggplot2:::plot_theme(p))$size
  }
  expect_equal(strip_pt(11), 14)
  expect_equal(strip_pt(22), 28)
})

test_that("return visibility", {
  # change 2: previously returned invisibly
  expect_true(withVisible(theme_dartR())$visible)
})

test_that("plots build without warnings", {
  p <- ggplot2::ggplot(data.frame(a = 1:10, g = rep(c("x", "y"), 5)),
                       ggplot2::aes(a, fill = g)) +
    ggplot2::geom_histogram(bins = 5) +
    ggplot2::facet_wrap(~g) +
    ggplot2::labs(title = "t", caption = "c", tag = "A") +
    theme_dartR()
  expect_no_warning(ggplot2::ggplotGrob(p))
})
