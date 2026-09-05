# Characterization tests for gl2eigenstrat — snapshot of CURRENT behaviour
# (review baseline, function-review campaign). Bugs are captured as-is and
# annotated; do not treat every expectation here as intended behaviour.

test_that("gl2eigenstrat writes the three EIGENSTRAT files for testset2.gl", {
  od <- file.path(tempdir(), "es_base")
  dir.create(od, showWarnings = FALSE)

  ret <- gl2eigenstrat(testset2.gl, outfile = "es", outpath = od, verbose = 0)
  expect_null(ret)

  g <- readLines(file.path(od, "es.eigenstratgeno"))
  s <- readLines(file.path(od, "es.snp"))
  ind <- readLines(file.path(od, "es.ind"))

  expect_length(g, nLoc(testset2.gl))          # one line per SNP
  expect_equal(nchar(g[1]), nInd(testset2.gl)) # one character per individual
  expect_length(s, nLoc(testset2.gl))
  expect_length(ind, nInd(testset2.gl))
  # defaults: chrom = 1, genetic position = 0 (pos.cM signature default,
  # despite the roxygen claiming [default 1]), physical position = 1
  expect_equal(s[1], "100049687-12-C/T 1 0 1 C T")
  expect_equal(ind[1], "AA010915 U Case")
})

test_that("gl2eigenstrat counts REFERENCE-allele copies in the geno file", {
  # [approved F1] the geno value is now 2 - score: EIGENSTRAT defines it as
  # copies of the reference allele (.snp column 5), while the dartR score
  # counts the alternate allele (.snp column 6). Matrix ind 1:5 x locus 1 is
  # 2 2 NA 2 2 (homozygous alternate) so the file line now starts 00900.
  od <- file.path(tempdir(), "es_geno")
  dir.create(od, showWarnings = FALSE)
  invisible(gl2eigenstrat(testset2.gl, outfile = "es2", outpath = od,
                          verbose = 0))
  g <- readLines(file.path(od, "es2.eigenstratgeno"))
  m <- as.matrix(testset2.gl)
  expect_equal(unname(m[1:5, 1]), c(2, 2, NA, 2, 2))
  expect_equal(substr(g[1], 1, 5), "00900")
  # NA is coded 9 throughout
  expect_equal(sum(vapply(g, function(l)
    lengths(regmatches(l, gregexpr("9", l))), 1L)), sum(is.na(m)))
})

test_that("gl2eigenstrat stops when no chromosome label can be encoded", {
  # [approved F2/F3] factor fields are coerced via as.character (previously
  # factor LEVEL CODES were written) and illegal chromosome values are
  # handled; the platypus labels (e.g. 'NC_041731.1_chromosome_4') are all
  # un-encodable, so the run now stops instead of writing garbage codes.
  od <- file.path(tempdir(), "es_chr")
  dir.create(od, showWarnings = FALSE)
  sub <- platypus.gl[1:5, 1:20]
  expect_error(
    gl2eigenstrat(sub, outfile = "plat",
                  snp.pos = "ChromPos_Platypus_Chrom_NCBIv1",
                  snp.chr = "Chrom_Platypus_Chrom_NCBIv1",
                  outpath = od, verbose = 0),
    "no locus has a chromosome value")
})

test_that("gl2eigenstrat maps X/Y/MT/XY labels and removes illegal-chromosome loci", {
  # [approved F2] as.character coercion + documented mapping (X -> 23,
  # Y -> 24, MT -> 90, XY -> 91); [approved F3] loci with illegal values
  # (non-numeric label, or < 1 such as 0) are removed with a warning at
  # verbose >= 1, and the geno file shrinks in step with the .snp file.
  od <- file.path(tempdir(), "es_map")
  dir.create(od, showWarnings = FALSE)
  sub <- testset2.gl[1:5, 1:8]
  lm <- sub@other$loc.metrics
  lm$chrom_test <- factor(c("1", "X", "Y", "0", "MT", "scaffold_7", "XY", "2"))
  sub@other$loc.metrics <- lm

  wout <- capture.output(invisible(
    gl2eigenstrat(sub, outfile = "map", snp.chr = "chrom_test",
                  outpath = od, verbose = 1)))
  expect_true(any(grepl("2 loci with illegal chromosome values", wout)))

  s <- readLines(file.path(od, "map.snp"))
  g <- readLines(file.path(od, "map.eigenstratgeno"))
  expect_length(s, 6)  # loci with chrom '0' and 'scaffold_7' removed
  expect_length(g, 6)
  ch_written <- as.numeric(sapply(strsplit(s, " +"), `[`, 2))
  expect_equal(ch_written, c(1, 23, 24, 90, 91, 2))
  expect_equal(sapply(strsplit(s, " +"), `[`, 1),
               locNames(sub)[c(1, 2, 3, 5, 7, 8)])
})

test_that("gl2eigenstrat is silent at verbose = 0", {
  od <- file.path(tempdir(), "es_sil")
  dir.create(od, showWarnings = FALSE)
  out <- capture.output(
    invisible(gl2eigenstrat(testset2.gl[1:5, 1:10], outfile = "sil",
                            outpath = od, verbose = 0)))
  expect_length(out, 0)
})

test_that("gl2eigenstrat rejects SilicoDArT data", {
  # [approved F4] accept = "SNP" now stops presence/absence data, which
  # previously produced a malformed 4-column .snp file (no loc.all, so the
  # ref/var columns vanished in a positional layout).
  expect_error(gl2eigenstrat(testset2.gs[1:5, 1:10], outfile = "gs",
                             outpath = tempdir(), verbose = 0),
               "SNP")
})

test_that("gl2eigenstrat validates sex.code and phen.value lengths", {
  # [approved F6] cbind previously recycled short vectors silently down the
  # .ind file (e.g. a 2-value sex.code alternating M/F).
  od <- file.path(tempdir(), "es_len")
  dir.create(od, showWarnings = FALSE)
  expect_error(gl2eigenstrat(testset2.gl[1:4, 1:5], outfile = "len",
                             sex.code = c("male", "female"),
                             outpath = od, verbose = 0),
               "sex.code")
  expect_error(gl2eigenstrat(testset2.gl[1:4, 1:5], outfile = "len",
                             phen.value = c("Case", "Control"),
                             outpath = od, verbose = 0),
               "phen.value")
  # full-length vectors are accepted
  expect_no_error(gl2eigenstrat(testset2.gl[1:4, 1:5], outfile = "len",
                                sex.code = rep("male", 4),
                                phen.value = rep("Control", 4),
                                outpath = od, verbose = 0))
})
