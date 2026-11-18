# Characterization tests for gl2genepop — snapshot of CURRENT behaviour (bugs
# included) captured pre-review at commit ed99203. Not a statement of correct
# behaviour; accepted diffs must map to approved review findings.

test_that("gl2genepop writes title, locus line, Pop blocks and 2-digit genotypes", {
  small <- gl.filter.allna(testset2.gl[1:6, 1:8], verbose = 0)
  # small is 6 individuals x 7 loci, 6 populations of 1
  td <- file.path(tempdir(), "gpop_char")
  dir.create(td, showWarnings = FALSE)
  out <- capture.output(res <- gl2genepop(small, outfile = "small.gen",
                                          outpath = td, verbose = 0))
  expect_s3_class(res, "data.frame")
  expect_named(res, "lines")
  f <- file.path(td, "small.gen")
  expect_true(file.exists(f))
  lns <- readLines(f)
  # title + locus line + (Pop + 1 individual) x 6 pops
  expect_length(lns, 2 + 2 * 6)
  expect_identical(lns[1], "Genepop output. Loci: 7 Populations: 6")
  expect_identical(lns[2], paste(locNames(small), collapse = ","))
  expect_identical(lns[3], "Pop")
  # allele coding A=01, T=02, C=03, G=04; missing 0000; id is pop_ind + comma
  expect_identical(lns[10],
    "EmmacMDBForb_AA010915, 0202 0101 0202 0101 0000 0303 0402")
  # pops are ordered alphabetically by default
  expect_identical(lns[4], "EmmacBurnBara_AA011723, 0202 0101 0202 0101 0202 0303 0402")
})

test_that("gl2genepop 3-digit format pads alleles and uses 000000 for missing", {
  small <- gl.filter.allna(testset2.gl[1:6, 1:8], verbose = 0)
  td <- file.path(tempdir(), "gpop_char")
  dir.create(td, showWarnings = FALSE)
  capture.output(gl2genepop(small, outfile = "small3.gen", outpath = td,
                            output.format = "3_digits", verbose = 0))
  lns <- readLines(file.path(td, "small3.gen"))
  expect_identical(lns[10],
    "EmmacMDBForb_AA010915, 002002 001001 002002 001001 000000 003003 004002")
})

test_that("gl2genepop on an all-monomorphic object writes one genotype per locus", {
  # Updated for approved finding F1: the locus list is now built with
  # unique() rather than negative indexing with duplicated(), so an
  # all-monomorphic object (one allele column per locus in the genind) no
  # longer collapses to two "0000" fields per row — each row carries one
  # correctly coded genotype per claimed locus.
  g3 <- gl.filter.allna(testset2.gl[1:9, 1:6], verbose = 0)
  pop(g3) <- rep(c("A", "B", "C"), each = 3)
  td <- file.path(tempdir(), "gpop_char")
  dir.create(td, showWarnings = FALSE)
  capture.output(gl2genepop(g3, outfile = "mono.gen", outpath = td,
                            verbose = 0))
  lns <- readLines(file.path(td, "mono.gen"))
  expect_identical(lns[1], "Genepop output. Loci: 5 Populations: 3")
  expect_identical(lns[4], "A_AA010915, 0202 0101 0202 0101 0000")
  # every data row carries exactly nLoc genotype fields
  data_rows <- grep(",", lns[-(1:2)], value = TRUE)
  expect_true(all(lengths(strsplit(sub("^[^,]+, ", "", data_rows), " ")) ==
                    nLoc(g3)))
})

test_that("gl2genepop pop.order errors on unlisted, misspelled or duplicated populations", {
  # Updated for approved finding F2: pop.order is now validated against
  # popNames(x); unlisted or misspelled names raise a fatal error naming
  # the offending entries instead of silently dropping populations.
  g3 <- gl.filter.allna(testset2.gl[1:9, 1:6], verbose = 0)
  pop(g3) <- rep(c("A", "B", "C"), each = 3)
  td <- file.path(tempdir(), "gpop_char")
  dir.create(td, showWarnings = FALSE)
  expect_error(
    capture.output(gl2genepop(g3, outfile = "sub.gen", outpath = td,
                              pop.order = c("C", "A"), verbose = 0)),
    "Missing from pop.order: B"
  )
  expect_error(
    capture.output(gl2genepop(g3, outfile = "typo.gen", outpath = td,
                              pop.order = c("C", "A", "Bee"), verbose = 0)),
    "Not found in the genlight object: Bee"
  )
  # a complete, correctly spelled ordering still works
  capture.output(gl2genepop(g3, outfile = "ord.gen", outpath = td,
                            pop.order = c("C", "A", "B"), verbose = 0))
  lns <- readLines(file.path(td, "ord.gen"))
  expect_true(any(grepl("^B_", lns)))
  first_of <- sapply(c("C_", "A_", "B_"), function(p) min(grep(p, lns)))
  expect_true(all(diff(first_of) > 0))  # populations appear in pop.order order
})

test_that("gl2genepop rejects SilicoDArT with a proper error and is silent at verbose = 0", {
  gs <- testset2.gs[1:4, 1:6]
  gs@other$loc.metrics <- testset2.gs@other$loc.metrics[1:6, ]
  # Updated for approved finding F5 (VRB2): the rejection now raises a
  # condition carrying the message (previously cat(error(...)) + bare stop()).
  expect_error(capture.output(gl2genepop(gs, verbose = 0)),
               "Only SNPs \\(diploid\\) data")
  # Updated for approved finding F3 (VRB5): the save-path message is gated
  # at verbose >= 2 (and no longer passes "\n" through file.path)
  small <- gl.filter.allna(testset2.gl[1:6, 1:8], verbose = 0)
  td <- file.path(tempdir(), "gpop_char")
  dir.create(td, showWarnings = FALSE)
  out <- capture.output(gl2genepop(small, outfile = "v0.gen", outpath = td,
                                   verbose = 0))
  expect_false(any(grepl("The genepop file is saved as", out)))
  expect_length(out, 0)
})

test_that("gl2genepop fails fast with a clear message on a single-individual object", {
  # Updated for approved finding F6: previously died inside the genind
  # conversion with "X is not a matrix".
  one <- gl.filter.allna(testset2.gl[1, 1:8], verbose = 0)
  td <- file.path(tempdir(), "gpop_char")
  dir.create(td, showWarnings = FALSE)
  expect_error(
    capture.output(gl2genepop(one, outfile = "one.gen", outpath = td,
                              verbose = 0)),
    "at least two individuals"
  )
})
