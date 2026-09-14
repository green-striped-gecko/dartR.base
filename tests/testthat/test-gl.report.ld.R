# Characterization tests for gl.report.ld
# Baseline snapshotted before review (review-gl.report.ld).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

prep_platypus_25 <- function() {
  x <- platypus.gl
  x <- gl.filter.callrate(x, threshold = 1, verbose = 0)
  pl <- x[, 1:25]
  pl@other$loc.metrics <- x@other$loc.metrics[1:25, ]
  pl
}

test_that("gl.report.ld baseline: all pooled pairs, 9 statistics", {
  pl <- prep_platypus_25()
  td <- file.path(tempdir(), "ldtest-report")
  dir.create(td, showWarnings = FALSE)
  out <- capture.output(
    lr <- gl.report.ld(pl, save = FALSE, outpath = td, nchunks = 2,
                       verbose = 0)
  )
  lr <- as.data.frame(lr)
  expect_equal(dim(lr), c(300, 9))
  expect_equal(colnames(lr),
               c("loc1", "loc2", "D", "Dprime", "r", "R2", "n", "X2", "p"))
  # deterministic first pair (loci 1 and 2, pooled over all 81 individuals)
  expect_equal(lr$D[1], -0.0021503531, tolerance = 1e-6)
  expect_equal(lr$R2[1], 2.580358e-03, tolerance = 1e-6)
  expect_equal(lr$n[1], 81)
  # [approved diff F6] no chunk files with save = FALSE
  expect_length(list.files(td, pattern = "LD_chunks"), 0)
  # [approved diff F5] fully silent at verbose = 0 (internal gl2gi silenced)
  expect_length(out, 0)
})

test_that("gl.report.ld rejects a genind", {
  # [approved diff F3, docs] genind objects were documented as accepted but
  # have always been rejected by the datatype check; the docs now say
  # genlight only
  pl <- prep_platypus_25()
  gi <- gl2gi(pl, verbose = 0)
  td <- file.path(tempdir(), "ldtest-report")
  expect_error(
    gl.report.ld(gi, save = FALSE, outpath = td, verbose = 0),
    "genind"
  )
})

test_that("gl.report.ld chunk restart returns cached results at any verbosity", {
  # [approved diff F2] the all-chunks-done early return no longer sits
  # inside if (verbose >= 2) (the rerun used to crash with "subscript out
  # of bounds" at verbose < 2), and chunk discovery now searches outpath
  # rather than the working directory
  pl <- prep_platypus_25()
  td <- file.path(tempdir(), "ldtest-restart")
  dir.create(td, showWarnings = FALSE)
  out1 <- capture.output(
    l1 <- gl.report.ld(pl, save = TRUE, outpath = td, nchunks = 2,
                       chunkname = "chr", verbose = 0)
  )
  expect_equal(nrow(l1), 300)
  expect_true(length(list.files(td, pattern = "LD_chunks_chr")) > 0)
  # rerun finds the chunks in outpath (working directory differs) and
  # returns the cached pairs without recomputing, at quiet verbosity
  out2 <- capture.output(
    l2 <- gl.report.ld(pl, save = TRUE, outpath = td, nchunks = 2,
                       chunkname = "chr", verbose = 0)
  )
  expect_equal(nrow(l2), 300)
  out3 <- capture.output(
    l3 <- gl.report.ld(pl, save = TRUE, outpath = td, nchunks = 2,
                       chunkname = "chr", verbose = 2)
  )
  expect_equal(nrow(l3), 300)
})

test_that("gl.report.ld rejects SilicoDArT data", {
  # [approved diff F7] presence/absence data used to be admitted and
  # LD-like statistics computed from 0/1 calls without warning
  gs <- testset.gs[, 1:10]
  gs@other$loc.metrics <- testset.gs@other$loc.metrics[1:10, ]
  td <- file.path(tempdir(), "ldtest-report")
  expect_error(
    capture.output(
      gl.report.ld(gs, save = FALSE, outpath = td, verbose = 0)
    ),
    "SilicoDArT"
  )
})
