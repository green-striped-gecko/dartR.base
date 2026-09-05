# Characterization tests for gl2gapit — snapshot of CURRENT behaviour (bugs
# included) captured pre-review at commit ed99203. Not a statement of correct
# behaviour; accepted diffs must map to approved review findings.

test_that("gl2gapit returns a hapmap-style data.frame with a duplicated header row", {
  small <- gl.filter.allna(testset.gl[1:6, 1:8], verbose = 0)
  # small is 6 individuals x 7 loci
  td <- file.path(tempdir(), "gapit_char")
  dir.create(td, showWarnings = FALSE)
  out <- capture.output(res <- gl2gapit(small, outfile = "gap", outpath = td,
                                        verbose = 0))
  expect_s3_class(res, "data.frame")
  # 1 header row + 7 loci; 11 hapmap columns + 6 individuals
  expect_equal(dim(res), c(nLoc(small) + 1, 11 + nInd(small)))
  # row 1 repeats the column names (GAPIT header=FALSE convention)
  expect_identical(unname(unlist(res[1, ])), colnames(res))
  expect_identical(colnames(res)[1:11],
                   c("rs", "alleles", "chrom", "pos", "strand", "assembly",
                     "center", "protLSID", "assayLSID", "panel", "QCcode"))
  # assembly is NA (was the hardcoded string "Oilpalm"; approved finding F6)
  expect_true(is.na(res$assembly[2]))
  # genotype coding: hom-ref doubles first allele, missing is "00"
  expect_identical(res[2, "AA010915"], "TT")   # locus 1, dosage 2, alleles C/T
  expect_identical(res[2, "AA032760"], "00")   # locus 1, NA
  expect_identical(res[3, "AA010915"], "AA")   # locus 2, dosage 0, alleles A/G
})

test_that("gl2gapit fills empty chromosome/position slots with dummies, silently at verbose = 0", {
  small <- gl.filter.allna(testset.gl[1:6, 1:8], verbose = 0)
  td <- file.path(tempdir(), "gapit_char")
  dir.create(td, showWarnings = FALSE)
  out <- capture.output(res <- gl2gapit(small, outfile = "gap2", outpath = td,
                                        verbose = 0))
  # the two slot warnings are now gated at verbose >= 2 (approved finding F4)
  expect_false(any(grepl("Chromosome slot is empty", out)))
  expect_false(any(grepl("Position slot is empty", out)))
  expect_length(out, 0)
  # dummy values: chrom "1" for all loci, pos 1..nLoc
  expect_identical(unique(res$chrom[-1]), "1")
  expect_identical(res$pos[-1], as.character(seq_len(nLoc(small))))
})

test_that("gl2gapit writes the hapmap table to outpath/outfile", {
  # Updated for approved finding F2 (DOC5): the function now writes the
  # tab-delimited hapmap file its parameters and message always promised.
  small <- gl.filter.allna(testset.gl[1:6, 1:8], verbose = 0)
  td <- file.path(tempdir(), "gapit_none")
  dir.create(td, showWarnings = FALSE)
  capture.output(res <- gl2gapit(small, outfile = "gap3", outpath = td,
                                 verbose = 0))
  f <- file.path(td, "gap3")
  expect_true(file.exists(f))
  lns <- readLines(f)
  # one line per data.frame row (header row + one per locus)
  expect_length(lns, nrow(res))
  expect_identical(strsplit(lns[1], "\t")[[1]], colnames(res))
})

test_that("gl2gapit recodes chromosome names to stable alphabetical codes", {
  # Updated for approved finding F3 (DOC5): codes are now assigned by
  # alphabetical order of the distinct names (stable across subsetting),
  # not by factor level indices. For this fixture the codes coincide with
  # the old behaviour.
  pl <- platypus.gl[1:4, 1:6]
  pl@other$loc.metrics <- platypus.gl@other$loc.metrics[1:6, ]
  pl$chromosome <- factor(c("chrB", "chrA", "chrC", "chrA", "chrB", "chrA"))
  pl$position <- c(10L, 20L, 30L, 40L, 50L, 60L)
  capture.output(res <- gl2gapit(pl, verbose = 0))
  expect_identical(res$chrom[-1], c("2", "1", "3", "1", "2", "1"))
  expect_identical(res$pos[-1], c("10", "20", "30", "40", "50", "60"))
})

test_that("gl2gapit rejects SilicoDArT input with a clear datatype error", {
  # Updated for approved finding F5 (DAT7): accept = "SNP" now rejects
  # presence/absence data up front instead of crashing in the fill loop.
  gs <- testset.gs[1:4, 1:6]
  gs@other$loc.metrics <- testset.gs@other$loc.metrics[1:6, ]
  expect_error(capture.output(gl2gapit(gs, verbose = 0)),
               "inappropriate object passed to function")
})
