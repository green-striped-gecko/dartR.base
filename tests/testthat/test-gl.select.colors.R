# Characterization tests for gl.select.colors
# Baseline snapshotted before review (review-gl.select.colors).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("the internal-caller contract holds", {
  pdf(NULL); on.exit(dev.off())
  o <- capture.output(
    cc <- gl.select.colors(library = "brewer", palette = "Blues",
                           select = c(7, 5), verbose = 0)
  )
  expect_length(o, 0)
  expect_identical(unname(cc), c("#2171B5", "#6BAED6"))
})

test_that("default path: silent, visible 9-colour vector", {
  pdf(NULL); on.exit(dev.off())
  o <- capture.output(
    v <- withVisible(gl.select.colors(plot.display = FALSE, verbose = 0))
  )
  expect_length(o, 0)
  expect_true(v$visible)   # interactive product; visible is correct
  expect_length(v$value, 9)
})

test_that("x= genlight drives ncolors and validates select length", {
  pdf(NULL); on.exit(dev.off())
  out <- capture.output(
    c8 <- gl.select.colors(x = testset.gl, library = "baseR",
                           palette = "rainbow", plot.display = FALSE,
                           verbose = 0)
  )
  expect_length(c8, nPop(testset.gl))
  expect_error(suppressWarnings(capture.output(
    gl.select.colors(x = testset.gl, library = "baseR", palette = "rainbow",
                     select = c(1, 2), plot.display = FALSE, verbose = 0)
  )))
})

test_that("invalid library value", {
  # [approved diff F1] an unknown library now coerces to the default
  # (scales/hue_pal) with a gated warning, instead of returning base R's
  # colors() FUNCTION via lexical scoping.
  pdf(NULL); on.exit(dev.off())
  out <- capture.output(
    r <- gl.select.colors(library = "bogus", plot.display = FALSE,
                          verbose = 0)
  )
  expect_length(out, 0)
  expect_true(is.character(r))
  expect_length(r, 9)
})

test_that("brewer count guarantees", {
  # [approved diff F2] requests below 3 are honoured exactly (trimmed from
  # the 3-colour pull); above-max requests return the palette maximum with
  # a gated warning.
  pdf(NULL); on.exit(dev.off())
  c2 <- gl.select.colors(library = "brewer", palette = "Blues", ncolors = 2,
                         plot.display = FALSE, verbose = 0)
  expect_length(c2, 2)
  out <- capture.output(
    c3 <- gl.select.colors(library = "brewer", palette = "Blues",
                           ncolors = 12, plot.display = FALSE, verbose = 0))
  expect_length(out, 0)   # warning gated at verbose >= 1
  expect_length(c3, 9)
})

test_that("select out of bounds", {
  # [approved diff F3] out-of-bounds indices are now a fatal error.
  pdf(NULL); on.exit(dev.off())
  expect_error(
    gl.select.colors(library = "brewer", palette = "Blues",
                     select = c(1, 15), plot.display = FALSE, verbose = 0),
    "select indices"
  )
})

test_that("baseR palette='heat'", {
  # [approved diff F4] 'heat' now dispatches to heat.colors.
  pdf(NULL); on.exit(dev.off())
  o <- capture.output(
    c6 <- gl.select.colors(library = "baseR", palette = "heat", ncolors = 6,
                           plot.display = FALSE, verbose = 2))
  expect_identical(unname(c6), grDevices::heat.colors(6))
  expect_false(any(grepl("not in Base R", o)))
})

test_that("datatype banner with x= at verbose 0", {
  # [approved diff F5] utils.check.datatype now receives verbose;
  # verbose = 0 is fully silent with x supplied.
  pdf(NULL); on.exit(dev.off())
  out <- capture.output(
    cx <- gl.select.colors(x = testset.gl, library = "baseR",
                           palette = "rainbow", plot.display = FALSE,
                           verbose = 0)
  )
  expect_length(out, 0)
})
