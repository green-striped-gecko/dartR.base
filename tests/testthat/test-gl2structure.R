# Characterization tests for gl2structure
# Baseline snapshotted pre-review (dev at ed99203). Assertions tagged
# [bug baseline] capture current defective behaviour on purpose; flip them
# only against an approved finding in the Phase A report.

str_sub <- function() testset2.gl[1:8, 1:15]

test_that("two-row-per-individual STRUCTURE layout with -9 missing", {
  td <- tempdir()
  f <- file.path(td, "str_base.str")
  if (file.exists(f)) file.remove(f)
  o <- capture.output(r <- gl2structure(str_sub(), outfile = "str_base.str",
                                        outpath = td, verbose = 0))
  expect_null(r)
  ln <- readLines(f)
  expect_length(ln, 1 + 2 * nInd(str_sub()))
  # header row = tab-separated locus names
  expect_equal(strsplit(ln[1], "\t")[[1]], locNames(str_sub()))
  # each individual contributes two adjacent rows labelled with its name
  labs <- vapply(strsplit(ln[-1], "\t"), `[`, "", 1)
  expect_equal(labs, rep(indNames(str_sub()), each = 2))
  # genotype coding: dosage 0 -> 1/1, 1 -> 1/2, 2 -> 2/2, NA -> -9/-9
  g <- as.matrix(str_sub())
  r1 <- as.integer(strsplit(ln[2], "\t")[[1]][-1])
  r2 <- as.integer(strsplit(ln[3], "\t")[[1]][-1])
  expected <- rbind(
    ifelse(is.na(g[1, ]), -9, ifelse(g[1, ] == 2, 2, 1)),
    ifelse(is.na(g[1, ]), -9, ifelse(g[1, ] >= 1, 2, 1)))
  expect_equal(r1, unname(expected[1, ]))
  expect_equal(r2, unname(expected[2, ]))
  expect_true(any(r1 == -9))
})

test_that("mismatched ind.names or add.columns are fatal", {
  expect_error(
    capture.output(gl2structure(str_sub(), ind.names = c("a", "b"),
                                outpath = tempdir(), verbose = 0)),
    "do not match")
  expect_error(
    capture.output(gl2structure(str_sub(), add.columns = data.frame(z = 1:3),
                                outpath = tempdir(), verbose = 0)),
    "does not match")
})

test_that("re-running without marker names overwrites the existing file", {
  # F1 fix applied: append = export.marker.names, so a headerless re-run
  # truncates instead of doubling the file
  td <- tempdir()
  f <- file.path(td, "str_app.str")
  if (file.exists(f)) file.remove(f)
  capture.output(gl2structure(str_sub(), export.marker.names = FALSE,
                              outfile = "str_app.str", outpath = td, verbose = 0))
  n1 <- length(readLines(f))
  capture.output(gl2structure(str_sub(), export.marker.names = FALSE,
                              outfile = "str_app.str", outpath = td, verbose = 0))
  n2 <- length(readLines(f))
  expect_equal(n1, 2 * nInd(str_sub()))
  expect_equal(n2, n1)   # F1 fix applied
})

test_that("SilicoDArT is rejected explicitly", {
  # F2 fix applied: accept = 'SNP' (DAT7) makes presence/absence data a
  # fatal error instead of a fabricated diploid STRUCTURE file
  td <- tempdir()
  f <- file.path(td, "str_gs.str")
  if (file.exists(f)) file.remove(f)
  expect_error(
    capture.output(gl2structure(testset2.gs[1:4, 1:10], outfile = "str_gs.str",
                                outpath = td, verbose = 0)),
    "SilicoDArT")
  expect_false(file.exists(f))
})

test_that("duplicated locus names keep one genotype column per locus", {
  # F3 fix applied: genotype columns are assigned positionally, so a name
  # collision no longer collapses columns: 15 names, 15 genotype columns
  td <- tempdir()
  f <- file.path(td, "str_dup.str")
  if (file.exists(f)) file.remove(f)
  dup <- str_sub()
  locNames(dup) <- c(rep("dupA", 2), locNames(str_sub())[3:15])
  capture.output(gl2structure(dup, outfile = "str_dup.str", outpath = td,
                              verbose = 0))
  ln <- readLines(f)
  expect_length(strsplit(ln[1], "\t")[[1]], 15)
  expect_length(strsplit(ln[2], "\t")[[1]], 16)  # 1 label + 15 loci (F3 fix)
  # the first duplicate's genotypes are its own, not the second's
  g <- as.matrix(dup)
  r1 <- as.integer(strsplit(ln[2], "\t")[[1]][-1])
  expect_equal(r1[1], unname(ifelse(is.na(g[1, 1]), -9,
                                    ifelse(g[1, 1] == 2, 2, 1))))
})

test_that("NULL return is invisible", {
  # F4 fix applied: invisible(NULL) rather than return(NULL)
  o <- capture.output(
    vis <- withVisible(gl2structure(str_sub(), outfile = "str_vis.str",
                                    outpath = tempdir(), verbose = 0)))
  expect_false(vis$visible)
})
