# Characterization tests for utils.heatmap
# Baseline snapshotted before review (branch review-utils.heatmap, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.

test_that("clustering and reordering match gplots::heatmap.2", {
  skip_if_not_installed("gplots")
  pdf(NULL); on.exit(dev.off())
  x <- as.matrix(mtcars)
  r1 <- utils.heatmap(x)
  r2 <- gplots::heatmap.2(x)
  expect_identical(r1$rowInd, r2$rowInd)
  expect_identical(r1$colInd, r2$colInd)
  expect_equal(dim(r1$carpet), dim(r2$carpet))
})

test_that("a matrix without dimnames plots", {
  pdf(NULL); on.exit(dev.off())
  # [approved diff H1] baseline: the added dendrogram auto-margin took
  # max(nchar(colnames(x))) on NULL names, sending -Inf into par(mar=)
  # ("invalid value specified for graphical parameter"); the gplots
  # original handles unnamed matrices.
  expect_error(utils.heatmap(matrix(rnorm(100), 10, 10)),
    NA)  # [approved diff H1]
})

test_that("fork features work: colored dendrogram labels and side colors", {
  pdf(NULL); on.exit(dev.off())
  x <- as.matrix(mtcars)
  expect_error(utils.heatmap(x, colRow = rep(c("red", "blue"), 16),
                             colCol = rep("purple", 11)), NA)
  expect_error(utils.heatmap(x, ColSideColors = rep("red", ncol(x))), NA)
  expect_error(utils.heatmap(x, ColSideColors = "red"), "ColSideColors")
})

test_that("scale='row' honours the return contract", {
  pdf(NULL); on.exit(dev.off())
  x <- as.matrix(mtcars)
  r <- utils.heatmap(x, scale = "row")
  expect_false(is.null(r$rowMeans))
  expect_false(is.null(r$rowSDs))
  expect_equal(dim(r$carpet), c(ncol(x), nrow(x)))
})

test_that("no label colours: no dendextend warning, plain labels", {
  # regression: without colRow/colCol (e.g. gl.plot.heatmap(D) with x = NULL)
  # dendextend::set("labels_col", NULL) warned "'x' is NULL so the result
  # will be NULL" for each dendrogram drawn
  pdf(NULL); on.exit(dev.off())
  x <- as.matrix(mtcars)
  expect_no_warning(utils.heatmap(x))
  expect_no_warning(utils.heatmap(x, dendrogram = "row"))
  expect_no_warning(utils.heatmap(x, dendrogram = "column"))
  # with label colours the dendrograms are still coloured without warnings
  expect_no_warning(utils.heatmap(x, colRow = rep(c("red", "blue"), 16),
                                  colCol = rep("purple", 11)))
})

test_that("gl.plot.heatmap without x gives no warning", {
  pdf(NULL); on.exit(dev.off())
  D <- dist(as.matrix(testset.gl[1:12, 1:50]))
  expect_no_warning(utils::capture.output(gl.plot.heatmap(D, verbose = 0)))
})
