# Characterization tests for utils.read.ped. Captured pre-review at
# upstream/dev ddaed27, then updated for the approved Phase C fixes
# (findings F1-F5, approved 2026-09-05); each changed expectation is
# annotated with its finding id. Findings and approval record:
# function-review/reports/dartR.base/utils.read.ped.md.
#
# utils.read.ped is a vendored copy of snpStats::read.pedfile; its only
# in-package caller is gl.report.ld.map.

toy_ped <- c(
  "FAM1 ind1 0 0 1 1 A A G G C C T T",
  "FAM1 ind2 0 0 2 1 A G G G C A 0 0",
  "FAM2 ind3 0 0 1 2 G G G T C A T T",
  "FAM2 ind4 0 0 2 2 A A T T C A T C")
toy_map <- c(
  "1 snp1 0 100",
  "1 snp2 0 200",
  "2 snp3 0 150",
  "2 snp4 0 300")

write_ped_fixture <- function(dir) {
  ped <- file.path(dir, "toy.ped")
  map <- file.path(dir, "toy.map")
  writeLines(toy_ped, ped)
  writeLines(toy_map, map)
  list(ped = ped, map = map)
}

ped_tmpdir <- function() {
  d <- tempfile("pedfix")
  dir.create(d)
  d
}

test_that("utils.read.ped parses a ped/map pair into genotypes, fam and map", {
  skip_if_not_installed("snpStats")
  d <- ped_tmpdir()
  f <- write_ped_fixture(d)

  r <- utils.read.ped(file = f$ped, snps = f$map)

  expect_named(r, c("genotypes", "fam", "map"))
  expect_s4_class(r$genotypes, "SnpMatrix")
  num <- as(r$genotypes, "numeric")
  # Numeric coding counts allele.2 = the SECOND allele encountered in file
  # order (not minor/major, not lexicographic).
  # snp1: a1 = A (first seen), a2 = G -> AA=0, AG=1, GG=2, AA=0
  expect_equal(unname(num[, "snp1"]), c(0, 1, 2, 0))
  # snp4: "0 0" honoured as missing under the default na.strings = "0"
  expect_equal(unname(num[, "snp4"]), c(0, NA, 0, 1))
  # fam columns
  expect_equal(r$fam$member, c("ind1", "ind2", "ind3", "ind4"))
  expect_equal(r$fam$sex, c(1, 2, 1, 2))
  expect_equal(r$fam$affected, c(1, 1, 2, 2))
  # map gains allele columns; the snps column is renamed
  expect_equal(r$map$allele.1, c("A", "G", "C", "T"))
  expect_equal(r$map$allele.2, c("G", "T", "A", "C"))
  expect_true("snp.names" %in% names(r$map))
})

test_that("utils.read.ped lex.order=TRUE swaps map alleles AND genotypes", {
  skip_if_not_installed("snpStats")
  d <- ped_tmpdir()
  f <- write_ped_fixture(d)

  rF <- utils.read.ped(file = f$ped, snps = f$map, lex.order = FALSE,
                       show_warnings = FALSE)
  rT <- utils.read.ped(file = f$ped, snps = f$map, lex.order = TRUE,
                       show_warnings = FALSE)

  # snp3 (C,A) and snp4 (T,C) are out of lexicographic order, so their map
  # alleles are swapped ...
  expect_equal(rT$map$allele.1, c("A", "G", "A", "C"))
  expect_equal(rT$map$allele.2, c("G", "T", "C", "T"))
  # [approved F1] ... and the genotypes now switch with them. Previously
  # switch.alleles' return value was discarded, so the genotype matrix was
  # byte-identical to the lex.order = FALSE run and dosages at snp3/snp4
  # counted the wrong allele relative to the returned map.
  expect_false(identical(rF$genotypes@.Data, rT$genotypes@.Data))
  numF <- as(rF$genotypes, "numeric")
  numT <- as(rT$genotypes, "numeric")
  # Unswapped loci are unchanged.
  expect_equal(numT[, "snp1"], numF[, "snp1"])
  expect_equal(numT[, "snp2"], numF[, "snp2"])
  # Swapped loci are complemented (count of the new allele.2):
  # snp3 was CC,CA,CA,CA -> 0,1,1,1 counting A; now 2,1,1,1 counting C.
  expect_equal(unname(numT[, "snp3"]), c(2, 1, 1, 1))
  # snp4 was TT,NA,TT,TC -> 0,NA,0,1 counting C; now 2,NA,2,1 counting T.
  expect_equal(unname(numT[, "snp4"]), c(2, NA, 2, 1))
})

test_that("utils.read.ped flags het-with-novel-allele as multi-allelic and NAs the locus", {
  skip_if_not_installed("snpStats")
  d <- ped_tmpdir()
  mp <- file.path(d, "multi.ped")
  # locus 1 alleles: A (a1, from i1), G (a2, from i2); i3 carries A/T where T
  # is a third allele
  writeLines(c("F1 i1 0 0 1 1 A A C C",
               "F1 i2 0 0 1 1 A G C C",
               "F1 i3 0 0 1 1 A T C A"), mp)

  # [approved F3] Previously the a1 match on the first allele masked the
  # novel second allele: A/T was coded 0 (homozygous A/A), the locus was
  # not flagged and no warning fired. Detection is now per allele column.
  expect_warning(r <- utils.read.ped(file = mp, show_warnings = TRUE),
                 "multi-allelic")
  num <- as(r$genotypes, "numeric")
  expect_true(all(is.na(num[, 1])))
  # The clean second locus is untouched.
  expect_equal(unname(num[, 2]), c(0, 0, 1))
})

test_that("utils.read.ped multi-allelic NA reset runs regardless of show_warnings", {
  skip_if_not_installed("snpStats")
  d <- ped_tmpdir()
  mp <- file.path(d, "multi3.ped")
  # i3 is T/T at locus 1: both alleles novel -> mallelic IS flagged
  writeLines(c("F1 i1 0 0 1 1 A A C C",
               "F1 i2 0 0 1 1 A G C C",
               "F1 i3 0 0 1 1 T T C A"), mp)

  expect_warning(rT <- utils.read.ped(file = mp, show_warnings = TRUE),
                 "multi-allelic")
  # show_warnings = FALSE stays quiet ...
  expect_no_warning(rF <- utils.read.ped(file = mp, show_warnings = FALSE))

  numT <- as(rT$genotypes, "numeric")
  numF <- as(rF$genotypes, "numeric")
  # Whole locus set to NA (documented snpStats behaviour).
  expect_true(all(is.na(numT[, 1])))
  # [approved F2] ... but the data cleaning now runs regardless.
  # Previously show_warnings = FALSE skipped the reset and i1/i2 kept
  # genotypes (0, 1, NA) at a locus classified as unreliable.
  expect_true(all(is.na(numF[, 1])))
})

test_that("utils.read.ped honours the split argument on data lines", {
  skip_if_not_installed("snpStats")
  d <- ped_tmpdir()
  cp <- file.path(d, "comma.ped")
  writeLines(c("F1,i1,0,0,1,1,A,A,C,C",
               "F1,i2,0,0,1,1,A,G,C,A"), cp)

  # [approved F4] Previously only the locus-count probe used split; the
  # per-line loop hardcoded "\t| +", so each whole line became field 1,
  # yielding line-string ids and all-NA genotypes with no error.
  r <- utils.read.ped(file = cp, split = ",")
  expect_equal(ncol(r$genotypes), 2)
  expect_equal(r$fam$member, c("i1", "i2"))
  num <- as(r$genotypes, "numeric")
  expect_equal(unname(num[, 1]), c(0, 1))
  expect_equal(unname(num[, 2]), c(0, 1))
})
