# Characterization tests for gl.read.fasta
# Baseline snapshotted before review (review-gl.read.fasta, Phase A).
# Phase C (2026-09-05): [BUG] pins for the approved findings flipped to the
# fixed behaviour, each tagged [approved F<n>]. F2/F5 are approved but their
# fix lands in the utils.read.fasta engine (its own review branch); the
# engine-dependent expectations here (locus sets, silence at verbose = 0)
# are therefore derived from the engine at run time, so this file holds on
# the wrapper branch alone and after the engine branch merges.

write_fasta <- function(lines, name) {
  path <- file.path(tempdir(), name)
  writeLines(lines, path)
  path
}

clean <- write_fasta(c(
  ">ind1", "ACACCAA-GT",
  ">ind2", "ACACTAACGT",
  ">ind3", "ACGCYCCTGT",
  ">ind4", "ACRCCCACGT",
  ">ind5", "ACACTGCCGT",
  ">ind6", "ACNCCASTGT"), "clean.fas")

clean2 <- write_fasta(c(
  ">ind1", "TAGAGC",
  ">ind2", "TAGATC",
  ">ind3", "TCGAKC",
  ">ind4", "TMGAGC",
  ">ind7", "TAGATC"), "clean2.fas")

mono <- write_fasta(c(
  ">m1", "ACGT",
  ">m2", "ACGT"), "mono.fas")

wrapped <- write_fasta(c(
  ">ind1", "ACACC", "AA-GT",
  ">ind2", "ACACT", "AACGT",
  ">ind3", "ACGCY", "CCTGT"), "wrapped.fas")

# the engine's view of the clean file: locus set and genotypes depend on
# which utils.read.fasta version is loaded, so derive rather than hardcode
engine_read <- function(f) {
  out <- capture.output(
    g <- dartR.base:::utils.read.fasta(f, parallel = FALSE, n.cores = NULL,
                                       verbose = 0))
  list(g = g, output = out)
}

test_that("default verbosity returns the genotype data", {
  # [approved F1] the fbm block ran gl.gen2fbm unconditionally and at
  # verbose <= 2 wiped @fbm too, returning a 0-individual shell; the
  # conversion now runs only when fbm = TRUE
  eng <- engine_read(clean)$g
  g <- suppressWarnings(gl.read.fasta(clean, verbose = 2))
  expect_equal(nInd(g), 6)
  expect_length(g@gen, 6)
  expect_true(is.null(g@fbm))
  m <- as.matrix(g)
  rownames(m) <- indNames(g)
  expect_equal(locNames(g), locNames(eng))
  # hand-verified column (C/T SNP with Y het), stable across engine versions
  expect_equal(unname(m[, "clean_5"]), c(0, 2, 1, 0, 2, 0))
})

test_that("verbose = 3 returns a plain object when fbm = FALSE", {
  # [approved F1] previously this path returned an FBM-backed object
  # although fbm was declined
  eng <- engine_read(clean)$g
  g <- suppressWarnings(gl.read.fasta(clean, verbose = 3))
  expect_true(is.null(g@fbm))
  expect_length(g@gen, 6)
  m <- as.matrix(g)
  rownames(m) <- indNames(g)
  expect_equal(nInd(g), 6)
  expect_equal(unique(ploidy(g)), 2)
  expect_equal(locNames(g), locNames(eng))
  expect_equal(unname(m[, "clean_5"]), c(0, 2, 1, 0, 2, 0))
  expect_null(g@position)                  # post-PR#330 convention
  expect_null(g@chromosome)
  # two entries: the gl.read.fasta call, then the internal
  # gl.recalc.metrics call (history is written before recalc runs)
  expect_length(g@other$history, 2)
  expect_false(is.null(g@other$ind.metrics))
  expect_equal(nrow(g@other$loc.metrics), nLoc(g))
})

test_that("fbm = TRUE returns an FBM-backed object", {
  # [approved F1] previously the data were lost at default verbosity
  # even when fbm was requested
  skip_if_not_installed("bigsnpr")
  g <- suppressWarnings(gl.read.fasta(clean, fbm = TRUE, verbose = 2))
  expect_length(g@gen, 0)
  expect_false(is.null(g@fbm))
  expect_equal(nInd(g), 6)
})

test_that("two files merge on individual names with NA fill", {
  eng1 <- engine_read(clean)$g
  eng2 <- engine_read(clean2)$g
  g <- suppressWarnings(gl.read.fasta(c(clean, clean2), verbose = 3))
  expect_equal(nInd(g), 7)
  expect_equal(nLoc(g), nLoc(eng1) + nLoc(eng2))
  expect_equal(indNames(g),
               c("ind1", "ind2", "ind3", "ind4", "ind5", "ind6", "ind7"))
  expect_equal(locNames(g), c(locNames(eng1), locNames(eng2)))
  m <- as.matrix(g)
  rownames(m) <- indNames(g)
  # individuals absent from a file are NA at that file's loci
  expect_true(all(is.na(m[c("ind5", "ind6"), locNames(eng2)])))
  expect_true(all(is.na(m["ind7", locNames(eng1)])))
  # hand-verified clean2 column (A/C SNP with M het): ind1..ind4, ind7;
  # stable across engine versions
  expect_equal(unname(m[c("ind1", "ind2", "ind3", "ind4", "ind7"),
                        "clean2_2"]), c(0, 0, 2, 1, 0))
  expect_equal(unique(ploidy(g)), 2)
  expect_equal(nrow(g@other$loc.metrics), nLoc(g))
})

test_that("a single all-monomorphic file is a clear fatal error", {
  # [approved F3] previously died downstream with an opaque
  # "argument of length 0"
  expect_error(
    capture.output(suppressWarnings(gl.read.fasta(mono, verbose = 0))),
    "no SNPs found")
})

test_that("line-wrapped FASTA is rejected with an informative error", {
  # [approved F4] previously mis-grouped records and died with
  # "missing value where TRUE/FALSE needed" (or mis-parsed silently)
  expect_error(suppressWarnings(gl.read.fasta(wrapped, verbose = 0)),
               "single line")
})

test_that("verbose = 0 is fully silent once engine messages are gated", {
  # [approved F5] the fix (gating the engine's multiallelic-skip and
  # no-polymorphism messages) lands in the utils.read.fasta review
  # branch; skip while the loaded engine still prints ungated so this
  # test activates when both branches are merged
  eng_out <- engine_read(clean)$output
  skip_if(length(eng_out) > 0,
          "utils.read.fasta engine messages not yet gated (sibling PR)")
  o <- capture.output(g <- gl.read.fasta(clean, verbose = 0))
  expect_length(o, 0)
})
