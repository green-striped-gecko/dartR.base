# Characterization tests for gl2paup.parsimony
# Baseline snapshotted pre-review (dev at ed99203). Assertions capture
# CURRENT behaviour, defects included, so Phase C changes surface as diffs.

make_gs <- function() {
  gs <- testset.gs[1:20, 1:100]
  gs@other$loc.metrics <- gs@other$loc.metrics[1:100, ]
  gs
}

test_that("standard nexus output structure is stable", {
  gs <- make_gs()
  td <- tempdir()
  o <- capture.output(rv <- gl2paup.parsimony(gs, outfileprefix = "chpars",
        outpath = td, nreps = 100, nbootstraps = 1000, verbose = 0))
  f <- file.path(td, "chpars_bootstrap.nex")
  expect_true(file.exists(f))
  ln <- readLines(f)
  expect_length(ln, 49)
  expect_equal(ln[1], "#NEXUS")
  expect_equal(ln[3], "     dimensions ntax = 20 nchar = 100 ;")
  expect_equal(ln[4], "     format datatype = standard gap = - ;")
  # first data record: individual name, spacing, 100-character 0/1/? string
  expect_match(ln[7], "^AA019237 +[01?]{100}$")
  expect_true(any(grepl("taxpartition pops =", ln, fixed = TRUE)))
  expect_true(any(grepl("bootstrap nreps=1000", ln, fixed = TRUE)))
  expect_null(rv)
})

test_that("verbose = 0 is fully silent with an invisible NULL", {
  gs <- make_gs()
  td <- tempdir()
  # F7 (warnings gated) and F8 (invisible NULL): no output at verbose = 0
  o <- capture.output(gl2paup.parsimony(gs, outfileprefix = "chv0",
        outpath = td, nreps = 1, nbootstraps = 10, verbose = 0))
  expect_length(o, 0)
})

test_that("bash mode: divisibility failure stops; scripts land in outpath", {
  gs <- make_gs()
  td <- tempdir()
  cwd <- file.path(td, "pars-cwd")
  dir.create(cwd, showWarnings = FALSE)
  old <- setwd(cwd)
  on.exit(setwd(old), add = TRUE)
  # F1: 10 %% 3 != 0 is now a fatal stop(error(...))
  expect_error(capture.output(gl2paup.parsimony(gs, outfileprefix = "chbash",
        outpath = td, out.type = "bash", nbootstraps = 10, ncpus = 3,
        base.dir.name = "/g/data/test", nreps = 1, verbose = 0)),
        "must be a multiple")
  # A divisible run writes all three generator scripts to outpath, not the cwd (F2)
  capture.output(gl2paup.parsimony(gs, outfileprefix = "chbash",
        outpath = td, out.type = "bash", nbootstraps = 10, ncpus = 5,
        base.dir.name = "/g/data/test", nreps = 1, verbose = 0))
  expect_setequal(list.files(td, pattern = "^generator_chbash"),
    c("generator_chbash_bootstraps.sh", "generator_chbash_consensus.sh",
      "generator_chbash_maketrees.sh"))
  expect_length(list.files(cwd, pattern = "^generator_chbash"), 0)
  ln <- readLines(file.path(td, "chbash_bootstrap.nex"))
  expect_true(any(grepl("bootstrap nreps=2 ", ln, fixed = TRUE)))
  # F3: maketrees jobs cd into base.dir.name, not a hardcoded personal directory
  mk <- readLines(file.path(td, "generator_chbash_maketrees.sh"))
  expect_false(any(grepl("ag3760", mk, fixed = TRUE)))
  expect_true(any(grepl("cd /g/data/test", mk, fixed = TRUE)))
  # F3: storage directive lists the project of base.dir.name plus if89
  expect_true(any(grepl("#PBS -l storage=gdata/xl04+gdata/test+gdata/if89",
      mk, fixed = TRUE)))
  # F5: consensus script starts with a shebang and has a fixed job name
  cons <- readLines(file.path(td, "generator_chbash_consensus.sh"))
  expect_equal(cons[1], "#!/bin/bash")
  expect_equal(cons[2], "#PBS -P xl04")
  expect_true(any(grepl("#PBS -N chbash_consensus", cons, fixed = TRUE)))
  expect_false(any(grepl("job${i}", cons, fixed = TRUE)))
})

test_that("SNP data is rejected; a single population converts cleanly", {
  td <- tempdir()
  expect_error(capture.output(gl2paup.parsimony(testset.gl[1:5, 1:50],
      outpath = td, verbose = 0)), "SilicoDArT")
  gs <- make_gs()
  pop(gs) <- rep("onlypop", nInd(gs))
  # F4: single-population input no longer crashes on the taxpartition loop
  capture.output(gl2paup.parsimony(gs, outfileprefix = "ch1p",
      outpath = td, nreps = 1, nbootstraps = 10, verbose = 0))
  ln <- readLines(file.path(td, "ch1p_bootstrap.nex"))
  expect_true(any(grepl("onlypop : 1-20;", ln, fixed = TRUE)))
})
