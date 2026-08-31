# Characterization tests for gl.report.pa
# Baseline snapshotted before review (branch review-gl.report.pa, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.

test_that("pairwise table matches an independent recomputation", {
  pn <- popNames(testset.gl)
  o <- capture.output(r <- gl.report.pa(testset.gl, verbose = 0))
  expect_equal(dim(r), c(435, 15))
  p1 <- as.matrix(testset.gl[pop(testset.gl) == pn[1], ])
  p2 <- as.matrix(testset.gl[pop(testset.gl) == pn[2], ])
  a1 <- colMeans(p1, na.rm = TRUE) / 2
  a2 <- colMeans(p2, na.rm = TRUE) / 2
  expect_equal(r$priv1[1],
               sum((a2 == 0 & a1 != 0) | (a2 == 1 & a1 != 1), na.rm = TRUE))
  expect_equal(r$priv2[1],
               sum((a1 == 0 & a2 != 0) | (a1 == 1 & a2 != 1), na.rm = TRUE))
  expect_equal(r$fixed[1], sum(abs(a1 - a2) == 1, na.rm = TRUE))
})

test_that("verbose 0 is silent", {
  # [approved diff F1] baseline: 32 warning lines leaked from gl.keep.loc
  # inside utils.pa.Chao ("no loci listed to keep") for zero-private pairs.
  o <- capture.output(r <- gl.report.pa(testset.gl, verbose = 0))
  expect_equal(length(o), 0)  # [approved diff F1]
})

test_that("Chao estimates depend only on the pair", {
  # [approved diff F1] baseline: utils.pa.Chao computed maf across ALL
  # populations in the dataset, so Chao values for a pair changed when
  # unrelated populations were present (full-x: 0/0; pair-only: 0/4), and
  # zero-private pairs fell back to the full locus set.
  pn <- popNames(testset.gl)
  pair <- gl.keep.pop(testset.gl, pop.list = pn[1:2], verbose = 0)
  o1 <- capture.output(rfull <- gl.report.pa(testset.gl, verbose = 0))
  o2 <- capture.output(rpair <- gl.report.pa(pair, verbose = 0))
  expect_equal(rfull$Chao1[1], rpair$Chao1[1])  # [approved diff F1]
  expect_equal(rfull$Chao2[1], rpair$Chao2[1])
})

test_that("plot.file without plot.display works", {
  # [approved diff F2] baseline: crashed with "object 'p3' not found"
  # after completing the whole computation.
  expect_error(
    gl.report.pa(testset.gl, plot.file = "patest", plot.dir = tempdir(),
                 verbose = 0),
    NA)  # [approved diff F2]
})

test_that("a mistyped method fails informatively", {
  # [approved diff F3] baseline: "object 'pall' not found".
  expect_error(
    gl.report.pa(testset.gl, method = "One2Rest", verbose = 0),
    "pairwise")  # [approved diff F3]
})

test_that("one2rest runs and its Sankey works with the default palette", {
  o <- capture.output(r <- gl.report.pa(testset.gl, method = "one2rest",
                                        verbose = 0))
  expect_equal(nrow(r), nPop(testset.gl))
  # [approved diff F4] baseline: plot.display = TRUE with the default
  # palette crashed ("argument is of length zero") -- the NULL-palette
  # default was handled in the pairwise branch only.
  pdf(NULL); on.exit(dev.off())
  expect_error(
    capture.output(gl.report.pa(testset.gl, method = "one2rest",
                                plot.display = TRUE, verbose = 0)),
    NA)  # [approved diff F4]
})

test_that("SilicoDArT report runs", {
  o <- capture.output(r <- gl.report.gs <- gl.report.pa(testset.gs,
                                                        verbose = 0))
  expect_equal(nrow(r), choose(nPop(testset.gs), 2))
  expect_true(is.na(r$Chao1[1]))
})
