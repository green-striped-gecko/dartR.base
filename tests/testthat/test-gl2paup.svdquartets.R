# Characterization tests for gl2paup.svdquartets
# Baseline snapshotted pre-review (dev at ed99203). Assertions capture
# CURRENT behaviour, defects included, so Phase C changes surface as diffs.

make_gl <- function() {
  gg <- testset2.gl[1:20, 1:100]
  gg@other$loc.metrics <- gg@other$loc.metrics[1:100, ]
  gg
}

make_gs <- function() {
  gs <- testset2.gs[1:20, 1:100]
  gs@other$loc.metrics <- gs@other$loc.metrics[1:100, ]
  gs
}

test_that("method 2 nexus structure and ambiguity coding are stable", {
  gg <- make_gl()
  td <- tempdir()
  o <- capture.output(rv <- gl2paup.svdquartets(gg, outfile = "chsvd2.nex",
        outpath = td, nbootstraps = 100, verbose = 0))
  ln <- readLines(file.path(td, "chsvd2.nex"))
  expect_length(ln, 55)
  expect_equal(ln[3], "     dimensions ntax = 20 nchar = 100 ;")
  expect_equal(ln[4], "     format datatype = dna gap = - ;")
  expect_match(ln[7], "^AA019237 +[ACGTRYSWKM?]{100}$")
  # het spot checks from the baseline run: G/T het -> K, A/G het -> R
  seq1 <- sub("^\\S+ +", "", ln[7])
  expect_equal(substr(seq1, 8, 8), "K")
  expect_equal(substr(seq1, 11, 11), "R")
  # F6: PAUP block filenames are derived from outfile
  expect_true(any(grepl("log file=chsvd2.txt;", ln, fixed = TRUE)))
  expect_true(any(grepl("treeFile=chsvd2.tre;", ln, fixed = TRUE)))
  expect_true(any(grepl("savetrees file=chsvd2_boot.tre", ln, fixed = TRUE)))
  expect_null(rv)
})

test_that("method 1 writes two lines per individual with split alleles", {
  gg <- make_gl()
  td <- tempdir()
  capture.output(gl2paup.svdquartets(gg, outfile = "chsvd1.nex",
        outpath = td, method = 1, nbootstraps = 100, verbose = 0))
  ln <- readLines(file.path(td, "chsvd1.nex"))
  expect_equal(ln[3], "     dimensions ntax = 40 nchar = 100 ;")
  expect_match(ln[7], "^AA019237_1 ")
  expect_match(ln[8], "^AA019237_2 ")
  # the G/T het that codes K under method 2 splits to G / T here
  s1 <- sub("^\\S+ +", "", ln[7])
  s2 <- sub("^\\S+ +", "", ln[8])
  expect_equal(substr(s1, 8, 8), "G")
  expect_equal(substr(s2, 8, 8), "T")
})

test_that("SilicoDArT output is written with datatype = standard", {
  gs <- make_gs()
  td <- tempdir()
  capture.output(gl2paup.svdquartets(gs, outfile = "chsvdgs.nex",
        outpath = td, nbootstraps = 100, verbose = 0))
  ln <- readLines(file.path(td, "chsvdgs.nex"))
  # F1: 0/1 presence-absence characters under a standard format declaration
  expect_equal(ln[4], "     format datatype = standard symbols = \"01\" gap = - ;")
  expect_match(ln[7], "^\\S+ +[01?]+$")
})

test_that("verbose = 0 is fully silent", {
  gg <- make_gl()
  td <- tempdir()
  # F2 (gl.sort quiet, warning gated) and F8 (invisible NULL): no output
  o <- capture.output(invisible(gl2paup.svdquartets(gg, outfile = "chsvdv0.nex",
        outpath = td, nbootstraps = 100, verbose = 0)))
  expect_length(o, 0)
})

test_that("method 0 is coerced with a warning; a single population converts cleanly", {
  gg <- make_gl()
  td <- tempdir()
  # F5: method = 0 is no longer silently accepted; warned at verbose >= 1
  o <- capture.output(gl2paup.svdquartets(gg, outfile = "chsvdm0.nex",
        outpath = td, method = 0, nbootstraps = 100, verbose = 1))
  ln <- readLines(file.path(td, "chsvdm0.nex"))
  expect_equal(ln[3], "     dimensions ntax = 20 nchar = 100 ;")  # one line per ind
  expect_true(any(grepl("method must be", o)))
  g1 <- make_gl()
  pop(g1) <- rep("onlypop", nInd(g1))
  # F4: single-population input no longer crashes on the taxpartition loop
  capture.output(gl2paup.svdquartets(g1, outfile = "chsvd1p.nex",
        outpath = td, verbose = 0))
  l1 <- readLines(file.path(td, "chsvd1p.nex"))
  expect_true(any(grepl("onlypop : 1-20;", l1, fixed = TRUE)))
})

test_that("a plain genlight passes the monomorph check; gl.sort still requires ind.metrics", {
  gg <- make_gl()
  plain <- new("genlight", as.matrix(gg))
  pop(plain) <- pop(gg)
  # F3: the monomorph flag check no longer crashes flag-less objects; the
  # remaining failure is gl.sort's own ind.metrics dependency (out of scope)
  expect_error(capture.output(gl2paup.svdquartets(plain,
        outfile = "chsvdpl.nex", outpath = tempdir(), verbose = 0)),
        "invalid subscript type")
})
