# Characterization tests for utils.read.fasta (internal, not exported)
# Baseline snapshotted before review (review-utils.read.fasta, Phase A).
# Phase C (2026-09-05): [BUG] pins for the approved findings F1-F9 flipped
# to the fixed behaviour; each flipped expectation carries an
# "[approved F<n>]" comment. LOW findings (F10, F11) were deferred; the
# NULL return pinned for F10 is unchanged in value.

urf <- dartR.base:::utils.read.fasta

write_fasta <- function(lines, name) {
  path <- file.path(tempdir(), name)
  writeLines(lines, path)
  path
}

# 10-column alignment, 6 individuals, 2-line records. Design (1-based cols):
#  col 3: A/G SNP, ind4 = R (het), ind6 = N (missing)
#  col 5: C/T SNP, ind3 = Y (het)
#  col 6: true triallelic A/C/G as three homozygote classes (A x3, C x2, G x1)
#  col 7: triallelic where the minority genotype is the het S=C/G -> detected
#  col 8: C/T SNP, ind1 = '-' (gap)
#  all other columns monomorphic
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

lowercase <- write_fasta(c(
  ">ind1", "CATAGA",
  ">ind2", "CaTAGA",
  ">ind3", "CGTaGA",
  ">ind4", "CATaGA"), "lowercase.fas")

unequal <- write_fasta(c(
  ">ind1", "ACGTACGTAC",
  ">ind2", "ACGTACGT",
  ">ind3", "ACGAACGTAC"), "unequal.fas")

mono <- write_fasta(c(
  ">m1", "ACGT",
  ">m2", "ACGT"), "mono.fas")

dupnames <- write_fasta(c(
  ">ind1", "ACGTAC",
  ">ind1", "ACGTTC",
  ">ind2", "ACGTAC"), "dupnames.fas")

test_that("SNP calling on a clean alignment: hets, ref/alt, locus naming", {
  g <- urf(clean, parallel = FALSE, n.cores = NULL, verbose = 0)
  expect_s4_class(g, "genlight")
  expect_equal(nInd(g), 6)
  # [approved F2] both triallelic columns (6: three homozygote classes;
  # 7: het-minority) are now detected from the allele pool and skipped
  expect_equal(nLoc(g), 3)
  expect_equal(locNames(g), c("clean_3", "clean_5", "clean_8"))
  expect_equal(g@loc.all, c("A/G", "C/T", "C/T"))
  expect_equal(unique(ploidy(g)), 2)
  # alignment column is encoded only in the locus name; the genome
  # coordinate slots stay NULL (post-PR#330 convention)
  expect_null(g@position)
  expect_null(g@chromosome)
  m <- as.matrix(g)
  # hand-verified: clean_5 = C/T with Y het -> 0,2,1,0,2,0 (correct path)
  expect_equal(unname(m[, "clean_5"]), c(0, 2, 1, 0, 2, 0))
  # ambiguity het R at clean_3/ind4 -> 1 (correct)
  expect_equal(unname(m["ind4", "clean_3"]), 1)
})

test_that("missing data (N and -) is coded NA", {
  # [approved F1] genotypes are classified explicitly; codes outside the
  # recognized genotype set (N, -, V/H/D/B) return NA, not a fake het
  g <- urf(clean, parallel = FALSE, n.cores = NULL, verbose = 0)
  m <- as.matrix(g)
  expect_true(is.na(m["ind6", "clean_3"]))   # N -> NA
  expect_true(is.na(m["ind1", "clean_8"]))   # gap -> NA
})

test_that("triallelic column of three homozygote classes is skipped", {
  # [approved F2] multiallelism is detected from the allele pool, not
  # from the two extreme genotype classes; column 6 (A/A x3, C/C x2,
  # G/G x1) is now dropped instead of surviving as A/G with fake hets
  g <- urf(clean, parallel = FALSE, n.cores = NULL, verbose = 0)
  expect_false("clean_6" %in% locNames(g))
})

test_that("multiallelic-skip message is gated: silent at 0, prints at 1", {
  # [approved F8] the skip and no-polymorphism messages are gated at
  # verbose >= 1 (they affect what the output contains, VRB4)
  o0 <- capture.output(g0 <- urf(clean, parallel = FALSE, n.cores = NULL,
                                 verbose = 0))
  expect_length(o0, 0)
  o1 <- capture.output(g1 <- urf(clean, parallel = FALSE, n.cores = NULL,
                                 verbose = 1))
  expect_true(any(grepl("more than 2 alleles", o1)))
  expect_true(any(grepl("6 7", o1)))         # the skipped columns, by name
})

test_that("lowercase (softmasked) bases are normalised, not fake alleles", {
  # [approved F3] sequences are toupper()ed on read: the case-mixed
  # monomorphic column is no longer a locus, and the case-mixed A/G
  # column scores the G/G individual as hom alt
  g <- urf(lowercase, parallel = FALSE, n.cores = NULL, verbose = 0)
  expect_equal(locNames(g), "lowercase_2")
  expect_equal(g@loc.all, "A/G")
  m <- as.matrix(g)
  expect_equal(unname(m[, "lowercase_2"]), c(0, 0, 2, 0))
})

test_that("unequal-length sequences are a fatal error naming the records", {
  # [approved F4] previously base-R recycling fabricated 8 loci from a
  # ragged alignment with only a generic matrix warning
  expect_error(urf(unequal, parallel = FALSE, n.cores = NULL, verbose = 0),
               "do not all have the same length")
  expect_error(urf(unequal, parallel = FALSE, n.cores = NULL, verbose = 0),
               "ind2")
})

test_that("no-polymorphism input returns NULL, silent at verbose = 0", {
  # [approved F8] the warning is gated; the NULL return is unchanged
  # (F10, LOW, deferred: the return is now written explicitly but its
  # value was NULL before and after)
  o <- capture.output(
    r <- urf(mono, parallel = FALSE, n.cores = NULL, verbose = 0))
  expect_null(r)
  expect_length(o, 0)
  o1 <- capture.output(
    r1 <- urf(mono, parallel = FALSE, n.cores = NULL, verbose = 1))
  expect_true(any(grepl("No polymorphism", o1)))
})

test_that("defaults are real: a bare call works", {
  # [approved F6] parallel = FALSE and verbose = NULL (resolved via
  # gl.check.verbosity) replace the self-referential defaults that died
  # with "promise already under evaluation"
  o <- capture.output(g <- urf(clean))
  expect_s4_class(g, "genlight")
  expect_equal(nLoc(g), 3)
})

test_that("reference allele follows allele frequency, not the modal genotype", {
  # [approved F5] ref/alt are chosen by summed allele counts from the
  # expanded genotypes: column 2 is Y,Y,Y,T,T,C (allele T = 7/12,
  # C = 5/12) so T is now the reference; previously the modal genotype
  # class (the het) flipped polarity to C/T
  refinv <- write_fasta(c(
    ">a1", "GYG", ">a2", "GYG", ">a3", "GYG",
    ">a4", "GTG", ">a5", "GTG", ">a6", "GCG"), "refinv.fas")
  g <- urf(refinv, parallel = FALSE, n.cores = NULL, verbose = 0)
  expect_equal(g@loc.all, "T/C")
  m <- as.matrix(g)
  expect_equal(unname(m[, 1]), c(1, 1, 1, 0, 0, 2))  # T/T -> 0, C/C -> 2
})

test_that("parallel path works: n.cores = NULL resolves, Windows falls back", {
  # [approved F7] n.cores = NULL resolves to parallel::detectCores();
  # on Windows the read falls back to serial (gated warning at
  # verbose >= 2) instead of erroring
  serial <- urf(clean, parallel = FALSE, n.cores = NULL, verbose = 0)
  p1 <- urf(clean, parallel = TRUE, n.cores = NULL, verbose = 0)
  expect_equal(as.matrix(p1), as.matrix(serial))
  p2 <- urf(clean, parallel = TRUE, n.cores = 2, verbose = 0)
  expect_equal(as.matrix(p2), as.matrix(serial))
})

test_that("duplicate individual names are preserved as-is", {
  g <- urf(dupnames, parallel = FALSE, n.cores = NULL, verbose = 0)
  expect_equal(indNames(g), c("ind1", "ind1", "ind2"))
})

test_that("duplicate names are refused before a multi-file merge", {
  # [approved F9] merge() joins equal names many-to-one; merge_gl_fasta
  # now refuses duplicated record names per file before joining
  g_dup <- urf(dupnames, parallel = FALSE, n.cores = NULL, verbose = 0)
  g_two <- urf(clean2, parallel = FALSE, n.cores = NULL, verbose = 0)
  expect_error(
    dartR.base:::merge_gl_fasta(list(g_dup, g_two), parallel = FALSE,
                                verbose = 0),
    "duplicate individual name")
})
