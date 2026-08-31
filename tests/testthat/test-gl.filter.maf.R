# Characterization tests for gl.filter.maf
# Baseline snapshotted before review (review-gl.filter.maf).
# Assertions marked [approved diff] were flipped in Phase C.

gl4 <- gl.keep.pop(testset.gl,
                   pop.list = c("EmmacBrisWive", "EmmacBurdMist",
                                "EmmacClarJack", "EmmacRussEube"),
                   verbose = 0)

test_that("global path matches an independent recomputation", {
  pdf(NULL); on.exit(dev.off())
  o <- capture.output(f <- gl.filter.maf(gl4, threshold = 0.05,
        plot.display = FALSE, verbose = 0))
  expect_length(o, 0)
  x1 <- utils.recalc.maf(gl4, verbose = 0)
  mymaf <- x1@other$loc.metrics$maf
  expect_equal(nLoc(f), length(which(mymaf >= 0.05)))   # 17 on gl4
  expect_equal(nLoc(f), 17)
  expect_equal(nrow(f@other$loc.metrics), nLoc(f))
  # NA-maf (all-NA) loci are removed by the global path
  expect_equal(sum(is.na(mymaf)), 3)
})

test_that("MAC interpretation for threshold > 1", {
  pdf(NULL); on.exit(dev.off())
  invisible(capture.output(f <- gl.filter.maf(gl4, threshold = 5,
        plot.display = FALSE, verbose = 0)))
  expect_equal(nLoc(f), 15)
})

test_that("by.pop path counts and history", {
  pdf(NULL); on.exit(dev.off())
  xcopy <- gl4
  invisible(capture.output(f <- gl.filter.maf(gl4, threshold = 0.05,
        by.pop = TRUE, pop.limit = 1, ind.limit = 5,
        plot.display = FALSE, verbose = 0)))
  expect_equal(nLoc(gl4) - nLoc(f), 221)
  expect_identical(xcopy, gl4)
  h <- f@other$history
  expect_equal(deparse(h[[length(h)]][[1]]), "gl.filter.maf")
})

test_that("by.pop with no qualifying loci or populations", {
  pdf(NULL); on.exit(dev.off())
  # [approved diff F1] baseline: when no populations met ind.limit (or
  # no loci qualified), x[,-integer(0)] subset to ZERO loci and
  # crashed. Now the object is returned unchanged, silent at verbose 0.
  o <- capture.output(f <- gl.filter.maf(gl4, threshold = 0.05,
        by.pop = TRUE, ind.limit = 50, plot.display = FALSE,
        verbose = 0))
  expect_length(o, 0)
  expect_equal(nLoc(f), nLoc(gl4))
})

test_that("plot.file with display off", {
  pdf(NULL); on.exit(dev.off())
  # [approved diff F2] baseline: crashed with "object 'p3' not found".
  # Now a gated notice; nothing saved; silent at verbose 0.
  o <- capture.output(f <- gl.filter.maf(gl4, threshold = 0.05,
        plot.display = FALSE, plot.file = "maf-x", verbose = 0))
  expect_length(o, 0)
  expect_equal(nLoc(f), 17)
})

test_that("threshold coercion warning at verbose 0; visibility", {
  pdf(NULL); on.exit(dev.off())
  # [approved diff F4] baseline: the range warning printed at verbose 0.
  o <- capture.output(f <- gl.filter.maf(gl4, threshold = 0.9,
        plot.display = FALSE, verbose = 0))
  expect_length(o, 0)
  # [approved diff F5] baseline: visible return.
  v <- withVisible(gl.filter.maf(gl4, threshold = 0.05,
        plot.display = FALSE, verbose = 0))
  expect_false(v$visible)
})

test_that("SilicoDArT data rejected", {
  pdf(NULL); on.exit(dev.off())
  # [approved diff F6] baseline: SilicoDArT objects ran through
  # utils.recalc.maf producing meaningless MAF values; now fatal.
  expect_error(capture.output(gl.filter.maf(testset.gs, threshold = 0.05,
        plot.display = FALSE, verbose = 0)))
})
