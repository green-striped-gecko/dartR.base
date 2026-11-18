# Characterization tests for gl2related
# Baseline snapshotted pre-review (dev at ed99203). Assertions capture
# CURRENT behaviour, defects included, so Phase C changes surface as diffs.

make_sub <- function() {
  gg <- testset2.gl[1:10, 1:20]
  gg@other$loc.metrics <- gg@other$loc.metrics[1:20, ]
  gg
}

test_that("two-column allele coding is stable (alleles 1/3, missing 0)", {
  gg <- make_sub()
  td <- tempdir()
  o <- capture.output(gtd <- gl2related(gg, outfile = "chrel.txt",
                                        outpath = td, verbose = 0))
  expect_s3_class(gtd, "data.frame")
  expect_equal(dim(gtd), c(10, 41))  # id column + 2 alleles x 20 loci
  gm <- as.matrix(gg)
  pick <- function(score) which(gm == score, arr.ind = TRUE)[1, ]
  cols <- function(j) c(2 * j, 2 * j + 1)
  h <- pick(1)
  expect_equal(unlist(gtd[h[1], cols(h[2])], use.names = FALSE), c(1, 3))
  r0 <- pick(0)
  expect_equal(unlist(gtd[r0[1], cols(r0[2])], use.names = FALSE), c(1, 1))
  r2 <- pick(2)
  expect_equal(unlist(gtd[r2[1], cols(r2[2])], use.names = FALSE), c(3, 3))
  na <- which(is.na(gm), arr.ind = TRUE)[1, ]
  expect_equal(unlist(gtd[na[1], cols(na[2])], use.names = FALSE), c(0, 0))
})

test_that("output file layout: tab-separated, no header, unquoted names", {
  gg <- make_sub()
  td <- tempdir()
  capture.output(gl2related(gg, outfile = "chrel2.txt", outpath = td,
                            verbose = 0))
  ln <- readLines(file.path(td, "chrel2.txt"))
  expect_length(ln, 10)
  # F1: quote = FALSE writes the name column without quote characters
  expect_match(ln[1], "^AA010915\t")
  expect_length(strsplit(ln[1], "\t")[[1]], 41)
})

test_that("save = FALSE writes nothing and verbose = 0 is silent", {
  gg <- make_sub()
  td2 <- file.path(tempdir(), "relnosave")
  dir.create(td2, showWarnings = FALSE)
  o <- capture.output(gtd <- gl2related(gg, outfile = "none.txt",
        outpath = td2, save = FALSE, verbose = 0))
  expect_false(file.exists(file.path(td2, "none.txt")))
  expect_length(o, 0)
  expect_s3_class(gtd, "data.frame")
})

test_that("SilicoDArT input is rejected", {
  # F2: accept = "SNP" stops silico presence scores being silently written
  # as heterozygote pairs
  gs <- testset2.gs[1:10, 1:20]
  gs@other$loc.metrics <- gs@other$loc.metrics[1:20, ]
  expect_error(capture.output(gl2related(gs, save = FALSE, verbose = 0)),
               "expecting SNP")
})
