# Characterization tests for gl.plot.heatmap
# Baseline snapshotted before review (dev_luis at 250b840, origin/dev
# f9ee087 merged). Assertions tagged [approved diff, change n] were flipped
# in Phase C to reflect the approved behaviour changes (report:
# function-review/reports/dartR.base/gl.plot.heatmap.md).

gl12 <- testset.gl[1:12, ]
D12 <- dist(as.matrix(gl12))

quiet_plot <- function(...) {
  pdf(NULL)
  on.exit(dev.off())
  out <- capture.output(res <- suppressWarnings(gl.plot.heatmap(..., verbose = 0)))
  list(res = res, out = out)
}

test_that("dist input: returns the utils.heatmap list invisibly", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  v <- withVisible(suppressWarnings(gl.plot.heatmap(D12, verbose = 0)))
  expect_false(v$visible)
  res <- v$value
  expect_type(res, "list")
  expect_named(res, c("rowInd", "colInd", "call", "carpet", "rowDendrogram",
                      "colDendrogram", "breaks", "col", "colorTable",
                      "layout"))
  expect_equal(dim(res$carpet), c(12, 12))
  # hclust on the same dist: reordering is deterministic
  expect_equal(res$rowInd, c(6, 1, 3, 4, 2, 11, 12, 5, 9, 10, 7, 8))
  expect_equal(res$colInd, res$rowInd)
  expect_equal(range(res$carpet), c(0, 4.873397), tolerance = 1e-6)
  expect_length(res$col, 255)
})

test_that("plot.out = FALSE draws nothing and returns NULL", {
  r <- quiet_plot(D12, plot.out = FALSE)
  expect_null(r$res)
})

test_that("input D is not modified", {
  D0 <- D12
  invisible(quiet_plot(D12, x = gl12, diag.na = TRUE))
  expect_identical(D12, D0)
})

test_that("verbose = 0 is silent", {
  # [approved diff, change 2] baseline: the default palette.divergent =
  # gl.colors("div") was evaluated at the session verbosity and printed
  # "Starting gl.colors", "Selected color type div", "Completed: gl.colors".
  r <- quiet_plot(D12)
  expect_length(r$out, 0)                                # [approved diff, change 2]
  r2 <- quiet_plot(D12, palette.divergent = gl.colors("div", verbose = 0))
  expect_length(r2$out, 0)
  r3 <- quiet_plot(D12, x = gl12)
  expect_length(r3$out, 0)
})

test_that("matrix input is plotted as supplied", {
  # [approved diff, change 1] baseline: as.dist() kept the lower triangle
  # only and zeroed the diagonal, so self-values and asymmetric cells
  # were discarded.
  m <- as.matrix(D12)
  diag(m) <- 5
  m[1, 2] <- 99
  r <- quiet_plot(m, dendrogram = "none")
  # the asymmetric cell reorders rows and columns differently: index the
  # plotted matrix back by name before comparing
  nm <- colnames(m)
  cc <- r$res$carpet[nm, nm]
  expect_equal(unname(diag(cc)), rep(5, 12))             # [approved diff, change 1]
  expect_true(any(cc == 99))                             # [approved diff, change 1]
  expect_equal(sum(cc == 99), 1)
  expect_true(all(cc == m) || all(cc == t(m)))
  # a matrix with the upper triangle empty is mirrored from the lower
  m2 <- as.matrix(D12)
  m2[upper.tri(m2)] <- NA
  r2 <- quiet_plot(m2, dendrogram = "none")
  expect_false(anyNA(r2$res$carpet))
  expect_equal(r2$res$carpet[rownames(m2), colnames(m2)], as.matrix(D12),
               ignore_attr = TRUE)
  # diag.na blanks the diagonal of the plotted matrix
  r3 <- quiet_plot(m, diag.na = TRUE, dendrogram = "none")
  expect_true(all(is.na(diag(r3$res$carpet[nm, nm]))))
  expect_equal(sum(is.na(r3$res$carpet)), 12)
})

test_that("non-square or name-mismatched matrix stops", {
  # [approved diff, change 1] baseline: as.dist() warned "non-square
  # matrix" and the function carried on with a truncated object.
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  m <- matrix(runif(20), 4, 5, dimnames = list(letters[1:4], letters[1:5]))
  expect_error(capture.output(gl.plot.heatmap(m, verbose = 0)),
               "square matrix")                          # [approved diff, change 1]
  m2 <- as.matrix(D12)
  rownames(m2) <- paste0("r", 1:12)
  expect_error(capture.output(gl.plot.heatmap(m2, verbose = 0)),
               "Row and column names")
})

test_that("fd input plots; x is ignored for fd", {
  fd <- gl.fixed.diff(gl12, verbose = 0)
  r <- quiet_plot(fd)
  expect_type(r$res, "list")
  expect_equal(dim(r$res$carpet), c(nPop(gl12), nPop(gl12)))
  # [approved diff, change 4] baseline: the population-colour block ran
  # as.matrix() on the fd list and stopped with "differing number of rows".
  r2 <- quiet_plot(fd, x = gl12)
  expect_type(r2$res, "list")                           # [approved diff, change 4]
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  out <- capture.output(invisible(gl.plot.heatmap(fd, x = gl12, verbose = 2)))
  expect_true(any(grepl("population-level matrix; x is ignored", out)))
})

test_that("legend labels and swatches pair, one row per population", {
  # [approved diff, change 3] baseline: legend colours were unique(colour),
  # so two populations sharing a colour left the last population with the
  # first colour (orange) instead of the one in its side bar (red).
  pal <- c("red", "red", "blue", "green", "orange", "purple", "black")
  seen <- NULL
  local_mocked_bindings(
    legend = function(x, y, legend, fill, ...) {
      seen <<- list(text = legend, fill = fill)
      invisible(NULL)
    },
    .package = "dartR.base")
  side <- NULL
  local_mocked_bindings(
    utils.heatmap = function(x, ColSideColors = NULL, ...) {
      side <<- ColSideColors
      list(carpet = x)
    },
    .package = "dartR.base")
  r <- quiet_plot(D12, x = gl12, palette_discrete = pal)
  expect_length(seen$text, 7)
  expect_length(seen$fill, 7)                           # [approved diff, change 3]
  expect_equal(seen$text, popNames(gl12))
  expect_equal(seen$fill, pal)                          # [approved diff, change 3]
  # side bar colours follow the individuals' populations in D's column order
  expect_equal(side,
               unname(setNames(pal, popNames(gl12))[as.character(pop(gl12))]))
  # a palette function is called with nPop(x)
  invisible(quiet_plot(D12, x = gl12, palette_discrete = rainbow))
  expect_equal(seen$fill, rainbow(7))
})

test_that("legend accepts a keyword position", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_error(capture.output(gl.plot.heatmap(D12, x = gl12,
                                              legendx = "topleft",
                                              verbose = 0)), NA)
})

test_that("par(mar) is restored after the legend", {
  # [approved diff, change 6] baseline: par(mar = c(1, 1, 1, 1)) before
  # legend() was never restored.
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  before <- par("mar")
  capture.output(invisible(gl.plot.heatmap(D12, x = gl12, verbose = 0)))
  expect_equal(par("mar"), before)                      # [approved diff, change 6]
})

test_that("x whose individuals do not match D: warning, no colours", {
  # [approved diff, change 5] baseline: fewer individuals in D than in x
  # dropped the colours with no message; same count with other names
  # stopped inside utils.heatmap ("ColSideColors must be ...").
  D8 <- dist(as.matrix(gl12[1:8, ]))
  seen <- NULL
  local_mocked_bindings(
    legend = function(x, y, legend, fill, ...) {
      seen <<- list(text = legend, fill = fill)
      invisible(NULL)
    },
    .package = "dartR.base")
  side <- "unset"
  local_mocked_bindings(
    utils.heatmap = function(x, ColSideColors = NULL, ...) {
      side <<- ColSideColors
      list(carpet = x)
    },
    .package = "dartR.base")
  # a subset of x's individuals is still coloured (every column matches)
  r <- quiet_plot(D8, x = gl12)
  expect_type(r$res, "list")
  expect_length(side, 8)                                # [approved diff, change 5]
  expect_equal(seen$text, popNames(gl12)[popNames(gl12) %in%
                                           as.character(pop(gl12))[1:8]])
  # renamed individuals: warning at verbose 1, plot drawn without colours
  m <- as.matrix(D12)
  dimnames(m) <- list(paste0("s", 1:12), paste0("s", 1:12))
  seen <- NULL
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  out <- capture.output(res <- gl.plot.heatmap(m, x = gl12, verbose = 1))
  expect_type(res, "list")                              # [approved diff, change 5]
  expect_true(any(grepl("12 of 12 columns of D are not individuals of x", out)))
  expect_null(seen)
  expect_null(side)
  out0 <- capture.output(invisible(gl.plot.heatmap(m, x = gl12, verbose = 0)))
  expect_length(out0, 0)
})

test_that("palette_discrete of the wrong length stops with a message", {
  # [approved diff, change 7] baseline: failed on names<- with "'names'
  # attribute [7] must be the same length as the vector [2]".
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_error(
    capture.output(gl.plot.heatmap(D12, x = gl12,
                                   palette_discrete = c("red", "blue"),
                                   verbose = 0)),
    "one per population of x")                          # [approved diff, change 7]
})

test_that("population-level dist with x: x ignored, plot still drawn", {
  Dp <- gl.dist.pop(testset.gl, verbose = 0)
  side <- "unset"
  local_mocked_bindings(
    utils.heatmap = function(x, ColSideColors = NULL, ...) {
      side <<- ColSideColors
      list(carpet = x)
    },
    .package = "dartR.base")
  r <- quiet_plot(Dp, x = testset.gl)
  expect_type(r$res, "list")
  expect_equal(dim(r$res$carpet), c(nPop(testset.gl), nPop(testset.gl)))
  expect_null(side)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  out <- capture.output(invisible(gl.plot.heatmap(Dp, x = testset.gl, verbose = 2)))
  expect_true(any(grepl("population-level matrix; x is ignored", out)))
  expect_false(any(grepl("Warning", out)))
})

test_that("verbose = 2 prints start, datatype and end lines", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  out <- capture.output(invisible(gl.plot.heatmap(D12, verbose = 2)))
  expect_equal(out[1], "Starting gl.plot.heatmap ")
  expect_equal(out[2], "  Processing a distance matrix")
  expect_equal(out[length(out)], "Completed: gl.plot.heatmap ")
  expect_length(out, 3)
})
