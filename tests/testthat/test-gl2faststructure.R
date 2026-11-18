# Characterization tests for gl2faststructure — snapshot of CURRENT behaviour
# (bugs included) captured pre-review at commit ed99203. Not a statement of
# correct behaviour; accepted diffs must map to approved review findings.

test_that("gl2faststructure writes 2 rows per individual, 6 dummy cols, 1/2/-9 coding", {
  small <- gl.filter.allna(testset2.gl[1:6, 1:8], verbose = 0)
  # small is 6 individuals x 7 loci
  td <- file.path(tempdir(), "fs_char")
  dir.create(td, showWarnings = FALSE)
  out <- capture.output(res <- gl2faststructure(small, outfile = "small.str",
                                                outpath = td, verbose = 0))
  expect_null(res)
  f <- file.path(td, "small.str")
  expect_true(file.exists(f))
  lns <- readLines(f)
  # two rows per diploid individual
  expect_length(lns, 2 * nInd(small))
  fields <- strsplit(lns, "\t")
  expect_true(all(lengths(fields) == nLoc(small) + 6))
  # first line snapshot (individual 1 = AA010915: 2 0 0 0 NA 0 1)
  expect_identical(lns[1], "1\t1\t1\t1\t1\t1\t2\t1\t1\t1\t-9\t1\t1")
  expect_identical(lns[2], "1\t1\t1\t1\t1\t1\t2\t1\t1\t1\t-9\t1\t2")
  # allele coding is 1/2 with -9 missing; dummy cols hold the row index, not
  # individual names
  geno <- unlist(lapply(fields, function(z) z[-(1:6)]))
  expect_true(all(geno %in% c("1", "2", "-9")))
  expect_identical(fields[[1]][1:6], rep("1", 6))
})

test_that("gl2faststructure het/hom recoding matches the 0/1/2 dosage input", {
  small <- gl.filter.allna(testset2.gl[1:6, 1:8], verbose = 0)
  td <- file.path(tempdir(), "fs_char")
  dir.create(td, showWarnings = FALSE)
  capture.output(gl2faststructure(small, outfile = "small2.str",
                                  outpath = td, verbose = 0))
  lns <- readLines(file.path(td, "small2.str"))
  m <- as.matrix(small)
  for (i in seq_len(nInd(small))) {
    r1 <- as.integer(strsplit(lns[2 * i - 1], "\t")[[1]][-(1:6)])
    r2 <- as.integer(strsplit(lns[2 * i], "\t")[[1]][-(1:6)])
    expected <- ifelse(is.na(m[i, ]), -18, ifelse(m[i, ] == 0, 2,
                       ifelse(m[i, ] == 2, 4, 3)))
    expect_equal(unname(r1 + r2), unname(expected))
  }
})

test_that("gl2faststructure rejects SilicoDArT input", {
  # Updated for approved finding F1 (DAT7): accept = "SNP" now rejects
  # presence/absence data instead of writing pseudo-diploid rows.
  gs <- testset2.gs[1:4, 1:6]
  gs@other$loc.metrics <- testset2.gs@other$loc.metrics[1:6, ]
  expect_error(
    capture.output(gl2faststructure(gs, outfile = "gs.str",
                                    outpath = tempdir(), verbose = 0))
  )
})

test_that("gl2faststructure returns NULL invisibly", {
  # Updated for approved finding F3 (FS10/VRB5): return(invisible(NULL)).
  small <- gl.filter.allna(testset2.gl[1:6, 1:8], verbose = 0)
  td <- file.path(tempdir(), "fs_char")
  dir.create(td, showWarnings = FALSE)
  vis <- withVisible(gl2faststructure(small, outfile = "vis.str",
                                      outpath = td, verbose = 0))
  expect_null(vis$value)
  expect_false(vis$visible)
})
