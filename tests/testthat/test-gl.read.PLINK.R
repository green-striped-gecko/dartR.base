# Characterization tests for gl.read.PLINK. Captured pre-review at
# upstream/dev ddaed27, then updated for the approved Phase C fixes
# (findings F1-F6, approved 2026-09-05); each changed expectation is
# annotated with its finding id. Findings and approval record:
# function-review/reports/dartR.base/gl.read.PLINK.md.
#
# Fixture: hand-built 4 ind x 4 SNP PLINK binary fileset. The .bed bytes were
# produced by PLINK 1.9 (v1.90b7.2) from the .ped/.map written below and are
# embedded verbatim so the test needs no plink executable.
#
# Truth table (alleles as written in the .ped):
#   snp1 (bim A1=G, A2=A): ind1 AA, ind2 AG, ind3 GG, ind4 AA
#   snp2 (bim A1=T, A2=G): ind1 GG, ind2 GG, ind3 GT, ind4 TT
#   snp3 (bim A1=A, A2=C): ind1 CC, ind2 CA, ind3 CA, ind4 CA
#   snp4 (bim A1=C, A2=T): ind1 TT, ind2 missing, ind3 TT, ind4 TC
# Dosage convention (unchanged; documented under F4): counts allele.2
# (bim col 6, the major allele), so hom-A2 = 2, hom-A1 = 0.

write_plink_fixture <- function(dir) {
  base <- file.path(dir, "toy")
  writeLines(c(
    "1\tsnp1\t0\t100\tG\tA",
    "1\tsnp2\t0\t200\tT\tG",
    "2\tsnp3\t0\t150\tA\tC",
    "2\tsnp4\t0\t300\tC\tT"), paste0(base, ".bim"))
  writeLines(c(
    "FAM1 ind1 0 0 1 1",
    "FAM1 ind2 0 0 2 1",
    "FAM2 ind3 0 0 1 2",
    "FAM2 ind4 0 0 2 2"), paste0(base, ".fam"))
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, 0xcb, 0x2f, 0xab, 0xb7)),
           paste0(base, ".bed"))
  base
}

plink_tmpdir <- function() {
  d <- tempfile("plinkfix")
  dir.create(d)
  d
}

expected_geno <- matrix(
  c(2, 2, 2, 2,
    1, 2, 1, NA,
    0, 1, 1, 2,
    2, 0, 1, 1),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("ind1", "ind2", "ind3", "ind4"),
                  c("snp1", "snp2", "snp3", "snp4")))

test_that("gl.read.PLINK at default-range verbosity returns genotypes in @gen", {
  skip_if_not_installed("snpStats")
  d <- plink_tmpdir()
  base <- write_plink_fixture(d)

  gl <- gl.read.PLINK(base, verbose = 0)

  expect_s4_class(gl, "dartR")
  expect_equal(indNames(gl), c("ind1", "ind2", "ind3", "ind4"))
  expect_equal(locNames(gl), c("snp1", "snp2", "snp3", "snp4"))
  # [approved F1] Previously length(gl@gen) == 0, nInd == 0 and
  # as.matrix() errored: gl.gen2fbm() always ran and the FBM copy was then
  # discarded at verbose <= 2. The object is now @gen-backed. (The fbm
  # slot is probed via the package's tolerant .fbm_or_null(): objects can
  # be built without the slot attribute.)
  expect_equal(nInd(gl), 4)
  expect_length(gl@gen, 4)
  expect_null(dartR.base:::.fbm_or_null(gl))
  expect_equal(as.matrix(gl), expected_geno)
})

test_that("gl.read.PLINK genotype content, slots and metrics at verbose = 3, fbm = FALSE", {
  skip_if_not_installed("snpStats")
  d <- plink_tmpdir()
  base <- write_plink_fixture(d)

  capture.output(gl <- gl.read.PLINK(base, verbose = 3, fbm = FALSE))
  # [approved F1] Previously verbose = 3 forced an FBM-backed object even
  # with fbm = FALSE (the flag was dead code); the backend now follows fbm.
  expect_null(dartR.base:::.fbm_or_null(gl))
  expect_gt(length(gl@gen), 0)

  m <- as.matrix(gl)
  # Dosage counts allele.2 (major); missing preserved as NA. Unchanged;
  # orientation now documented in @details (F4, documentation only).
  expect_equal(m, expected_geno)

  # Genome coordinates from the .bim land in the slots (PR #330 convention).
  expect_equal(as.integer(gl@position), c(100L, 200L, 150L, 300L))
  expect_equal(as.character(gl@chromosome), c("1", "1", "2", "2"))

  # .fam columns land in ind.metrics.
  expect_equal(as.character(gl@other$ind.metrics$Family),
               c("FAM1", "FAM1", "FAM2", "FAM2"))
  expect_equal(as.numeric(gl@other$ind.metrics$Sex), c(1, 2, 1, 2))

  # bim alleles land in loc.metrics; snp.name renamed AlleleID.
  expect_equal(as.character(gl@other$loc.metrics$AlleleID),
               c("snp1", "snp2", "snp3", "snp4"))
  expect_equal(as.character(gl@other$loc.metrics$allele.1),
               c("G", "T", "A", "C"))
  expect_equal(length(gl@other$history), 1L)
})

test_that("gl.read.PLINK fbm = TRUE returns a file-backed object", {
  skip_if_not_installed("snpStats")
  skip_if_not_installed("bigsnpr")
  d <- plink_tmpdir()
  base <- write_plink_fixture(d)

  # [approved F1] The fbm argument now decides the backend.
  gl <- gl.read.PLINK(base, fbm = TRUE, verbose = 0)
  expect_false(is.null(dartR.base:::.fbm_or_null(gl)))
  expect_equal(as.matrix(gl), expected_geno)
})

test_that("gl.read.PLINK ind.metafile rows are matched by id; pop applied", {
  skip_if_not_installed("snpStats")
  d <- plink_tmpdir()
  base <- write_plink_fixture(d)
  im <- file.path(d, "indmeta.csv")
  # ids deliberately NOT in .fam order
  write.csv(data.frame(id = c("ind3", "ind1", "ind4", "ind2"),
                       pop = c("P3", "P1", "P4", "P2"),
                       tag = c("t3", "t1", "t4", "t2")),
            im, row.names = FALSE)

  gl <- gl.read.PLINK(base, ind.metafile = im, verbose = 0)

  # [approved F2] Previously the rows were cbound in file order (ind1
  # received ind3's metadata) and the pop column was never applied to
  # pop(gl) (every individual stayed "A"). Rows are now matched by id and
  # pop is assigned.
  expect_equal(as.character(gl@other$ind.metrics$id[1]), "ind1")
  expect_equal(as.character(gl@other$ind.metrics$pop[1]), "P1")
  expect_equal(as.character(gl@other$ind.metrics$tag[1]), "t1")
  expect_equal(as.character(pop(gl)), c("P1", "P2", "P3", "P4"))
})

test_that("gl.read.PLINK loc.metafile without AlleleID column stops with the intended message", {
  skip_if_not_installed("snpStats")
  d <- plink_tmpdir()
  base <- write_plink_fixture(d)
  lm <- file.path(d, "locmeta.csv")
  write.csv(data.frame(marker = paste0("snp", 1:4), score = 1:4),
            lm, row.names = FALSE)

  # [approved F3] Previously the mandatory-column check printed but did not
  # stop; execution continued into an unrelated "undefined columns" error.
  expect_error(gl.read.PLINK(base, loc.metafile = lm, verbose = 0),
               "AlleleID column absent")
})

test_that("gl.read.PLINK merges a valid loc.metafile by AlleleID row name", {
  skip_if_not_installed("snpStats")
  d <- plink_tmpdir()
  base <- write_plink_fixture(d)
  lm <- file.path(d, "locmeta.csv")
  # rows deliberately scrambled; merge is by AlleleID so it must realign
  write.csv(data.frame(AlleleID = c("snp2", "snp1", "snp4", "snp3"),
                       myscore = c(20, 10, 40, 30)),
            lm, row.names = FALSE)

  gl <- gl.read.PLINK(base, loc.metafile = lm, verbose = 0)
  expect_equal(gl@other$loc.metrics$myscore, c(10, 20, 30, 40))
})

test_that("gl.read.PLINK is silent at verbose = 0 on the .bed path", {
  skip_if_not_installed("snpStats")
  d <- plink_tmpdir()
  base <- write_plink_fixture(d)
  out <- capture.output(gl <- gl.read.PLINK(base, verbose = 0))
  expect_length(out, 0)
})
