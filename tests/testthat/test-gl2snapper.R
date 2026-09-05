# Characterization tests for gl2snapper
# Baseline snapshotted pre-review (dev at ed99203). Assertions tagged
# [bug baseline] capture current defective behaviour on purpose; flip them
# only against an approved finding in the Phase A report.

snap_sub <- function() testset.gl[1:8, 1:15]

test_that("nexus output has the documented snapper structure", {
  td <- tempdir()
  o <- capture.output(r <- gl2snapper(snap_sub(), outfile = "snap_base.nex",
                                      outpath = td, verbose = 0))
  expect_length(o, 0)          # silent when no preprocessing options are used
  expect_null(r)
  ln <- readLines(file.path(td, "snap_base.nex"))
  expect_equal(ln[1], "#nexus")
  expect_equal(ln[2], "BEGIN DATA;")
  expect_match(ln[3], "ntax = 8 nchar = 15")
  expect_match(ln[4], 'integerdata missing=\\? symbols="012"')
  expect_equal(ln[5], "matrix")
  # one matrix row per individual, taxon label = pop_ind, genotypes 0/1/2/?
  rows <- ln[6:13]
  expect_length(rows, nInd(snap_sub()))
  expect_match(rows[1], "^EmmacMDBForb_AA010915  [012?]{15}$")
  # genotype string mirrors as.matrix with NA -> ?
  g <- as.matrix(snap_sub())[1, ]
  g[is.na(g)] <- "?"
  expect_equal(sub("^\\S+  ", "", rows[1]), paste(g, collapse = ""))
  expect_equal(ln[14], ";")
  expect_equal(ln[15], "end;")
})

test_that("SilicoDArT is rejected explicitly", {
  expect_error(
    capture.output(gl2snapper(testset.gs, outpath = tempdir(), verbose = 0)),
    "SilicoDArT")
})

test_that("individual names with spaces or duplicates are mangled to safe labels", {
  td <- tempdir()
  s <- testset.gl[1:6, 1:10]
  indNames(s)[2] <- indNames(s)[1]
  indNames(s)[3] <- "AA space"
  o <- capture.output(gl2snapper(s, outfile = "snap_names.nex",
                                 outpath = td, verbose = 0))
  # F1 fix applied: the name warnings are gated at verbose >= 2, silent here
  expect_length(o, 0)
  ln <- readLines(file.path(td, "snap_names.nex"))
  expect_true(any(grepl("_AA010915_1  ", ln)))
  expect_true(any(grepl("_AA010915_2  ", ln)))
  expect_true(any(grepl("_AA_space  ", ln)))
})

test_that("preprocessing options are silent at verbose 0", {
  td <- tempdir()
  set.seed(42)
  o <- capture.output(gl2snapper(testset.gl[, 1:100], outfile = "snap_pre.nex",
                                 outpath = td, rm.autapomorphies = TRUE,
                                 nloc = 20, verbose = 0))
  # F1 fix applied: verbose passed through to gl.allele.freq / gl.drop.loc
  # and the two cat(report()) calls gated at verbose >= 2
  expect_length(o, 0)
  ln <- readLines(file.path(td, "snap_pre.nex"))
  expect_match(ln[3], "nchar = 20 ")   # nloc subsample honoured
})

test_that("rm.autapomorphies default is FALSE", {
  # F2 fix applied: roxygen now documents [default FALSE], matching the
  # signature (signature unchanged)
  expect_false(formals(gl2snapper)$rm.autapomorphies)
})

test_that("NULL return is invisible", {
  # F5 fix applied: invisible(NULL) rather than return(NULL)
  td <- tempdir()
  o <- capture.output(
    vis <- withVisible(gl2snapper(snap_sub(), outfile = "snap_vis.nex",
                                  outpath = td, verbose = 0)))
  expect_false(vis$visible)
})
