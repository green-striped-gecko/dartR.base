# Characterization tests for gl2gds — snapshot of CURRENT behaviour (bugs
# included) captured pre-review at commit ed99203. Not a statement of correct
# behaviour; accepted diffs must map to approved review findings.

test_that("gl2gds default path writes a gds with zero pos/chrom and raw dosage", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  small <- gl.filter.allna(testset2.gl[1:6, 1:8], verbose = 0)
  # small is 6 individuals x 7 loci
  td <- file.path(tempdir(), "gds_char")
  dir.create(td, showWarnings = FALSE)
  f <- file.path(td, "small.gds")
  if (file.exists(f)) unlink(f)
  out <- capture.output(res <- gl2gds(small, outfile = "small.gds",
                                      outpath = td, verbose = 0))
  expect_null(res)
  expect_true(file.exists(f))
  gf <- SNPRelate::snpgdsOpen(f)
  on.exit(SNPRelate::snpgdsClose(gf), add = TRUE)
  # with default snp.chr = "0", snp.pos = "0" the sort is a no-op, so locus
  # order is preserved and pos/chrom are all zero
  expect_identical(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "snp.id")),
                   locNames(small))
  expect_identical(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "snp.position")),
                   rep(0L, nLoc(small)))
  expect_true(all(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "snp.chromosome")) == 0))
  expect_identical(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "snp.allele")),
                   unname(small@loc.all))
  # genotype is now stored as 2 - dosage (approved finding F2): SNPRelate
  # defines the stored value as the count of the FIRST allele of snp.allele,
  # while the genlight dosage counts the second; NA -> 3
  g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "genotype"))
  m <- t(as.matrix(small))
  m <- 2 - m
  m[is.na(m)] <- 3
  expect_equal(unname(g), unname(m))
  # with the orientation fixed, snpgdsSNPRateFreq agrees with the genlight
  # first-allele frequencies (approved finding F2)
  freq <- SNPRelate::snpgdsSNPRateFreq(gf)$AlleleFreq
  expected_freq <- 1 - colMeans(as.matrix(small), na.rm = TRUE) / 2
  expect_equal(unname(freq), unname(expected_freq))
})

test_that("gl2gds with snp.chr/snp.pos maps coordinates 1:1 to each snp.id", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")
  # Updated for approved finding F1 (BLOCKER): a single sort permutation is
  # now applied to genotype rows, snp.id, snp.allele, positions and
  # chromosomes, so every record carries its own coordinates (previously 7
  # of these 8 records carried another locus's position). Dosage inversion
  # per approved finding F2.
  pl <- platypus.gl[1:5, 1:8]
  pl@other$loc.metrics <- platypus.gl@other$loc.metrics[1:8, ]
  td <- file.path(tempdir(), "gds_char")
  dir.create(td, showWarnings = FALSE)
  f <- file.path(td, "perm.gds")
  if (file.exists(f)) unlink(f)
  capture.output(gl2gds(pl, outfile = "perm.gds", outpath = td,
                        snp.pos = "ChromPos_Platypus_Chrom_NCBIv1",
                        snp.chr = "Chrom_Platypus_Chrom_NCBIv1", verbose = 0))
  gf <- SNPRelate::snpgdsOpen(f)
  on.exit(SNPRelate::snpgdsClose(gf), add = TRUE)
  ids <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "snp.id"))
  pos <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "snp.position"))
  chrom <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "snp.chromosome"))
  alleles <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "snp.allele"))
  true_pos_all <- pl@other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
  true_chr_all <- as.character(
    pl@other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1)
  # records are sorted by (chromosome, position)
  expect_identical(ids,
                   locNames(pl)[order(true_chr_all, true_pos_all)])
  # every record carries its own position, chromosome and alleles
  idx <- match(ids, locNames(pl))
  expect_equal(pos, true_pos_all[idx])
  expect_identical(chrom, true_chr_all[idx])
  expect_identical(alleles, unname(pl@loc.all[idx]))
  # the genotype rows travel with snp.id, stored as 2 - dosage (F2)
  g <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "genotype"))
  m <- t(as.matrix(pl))
  m <- 2 - m
  m[is.na(m)] <- 3
  expect_equal(unname(g), unname(m[ids, ]))
})

test_that("gl2gds is silent at verbose = 0", {
  skip_if_not_installed("SNPRelate")
  # Updated for approved finding F3 (VRB5): the "Structure of gds file"
  # banner, snpgdsSummary and print(genofile) are now gated at verbose >= 3.
  small <- gl.filter.allna(testset2.gl[1:6, 1:8], verbose = 0)
  td <- file.path(tempdir(), "gds_char")
  dir.create(td, showWarnings = FALSE)
  f <- file.path(td, "v0.gds")
  if (file.exists(f)) unlink(f)
  out <- capture.output(gl2gds(small, outfile = "v0.gds", outpath = td,
                               verbose = 0))
  expect_false(any(grepl("Structure of gds file", out)))
  expect_length(out, 0)
})

test_that("gl2gds rejects SilicoDArT input", {
  skip_if_not_installed("SNPRelate")
  # Updated for approved finding F4 (DAT7): accept = "SNP" now rejects
  # presence/absence data instead of writing a pseudo-SNP gds.
  gs <- testset2.gs[1:4, 1:6]
  gs@other$loc.metrics <- testset2.gs@other$loc.metrics[1:6, ]
  td <- file.path(tempdir(), "gds_char")
  dir.create(td, showWarnings = FALSE)
  expect_error(
    capture.output(gl2gds(gs, outfile = "gs.gds", outpath = td, verbose = 0)),
    "inappropriate object passed to function"
  )
})
