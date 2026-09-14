# Characterization tests for gl2vcf
# Baseline snapshotted pre-review (dev at ed99203). Assertions tagged
# [bug baseline] capture current defective behaviour on purpose; flip them
# only against an approved finding in the Phase A report.
#
# gl2vcf shells out to PLINK 1.9. The coordinate-handling assertions below
# do not need a working PLINK: gl2plink writes the intermediate
# gl_plink_temp.map (into tempdir() since the F8 fix) before the plink
# call runs, so the POS values fed to plink are inspectable even when the
# binary cannot run (the call is wrapped in try()). A dummy plink file
# satisfies the up-front existence guard added under F6. The end-to-end
# VCF assertions are gated on the PLINK19_DIR environment variable naming
# a directory that contains plink[.exe].

vcf_fixture <- function() {
  sub <- dartR.data::platypus.gl[1:6, 1:8]
  sub@position <- NULL      # clear the stale tag-offset copy shipped in dartR.data
  sub@chromosome <- NULL
  sub
}

# a directory holding a dummy (non-runnable) plink, to get past the F6
# existence guard while still stopping short of a real PLINK run
dummy_plink_dir <- function(case) {
  td <- file.path(tempdir(), paste0("vcfmap_", case))
  dir.create(td, showWarnings = FALSE)
  file.create(file.path(td, "plink.exe"))
  td
}

run_gl2vcf_map <- function(x, case, ...) {
  td <- dummy_plink_dir(case)
  unlink(file.path(tempdir(), "gl_plink_temp.map"))
  try(suppressMessages(capture.output(
    gl2vcf(x, plink.bin.path = td,
           outfile = "out", outpath = td, verbose = 0, ...))),
    silent = TRUE)
  # F8 fix applied: the ped/map intermediates are written to tempdir(),
  # not outpath
  read.table(file.path(tempdir(), "gl_plink_temp.map"),
             col.names = c("chrom", "id", "cM", "pos"))
}

test_that("with a NULL position slot and no snp.pos, POS and CHROM are zero", {
  skip_if_not_installed("dartR.data")
  map <- run_gl2vcf_map(vcf_fixture(), "default")
  expect_equal(map$pos, rep(0L, 8))
  expect_equal(as.character(map$chrom), rep("0", 8))
  expect_equal(map$id, locNames(vcf_fixture()))
})

test_that("snp.pos/snp.chr fill the slots from loc.metrics when the slots are NULL", {
  skip_if_not_installed("dartR.data")
  x <- vcf_fixture()
  map <- run_gl2vcf_map(x, "named",
                        snp.pos = "ChromPos_Platypus_Chrom_NCBIv1",
                        snp.chr = "Chrom_Platypus_Chrom_NCBIv1")
  expect_equal(map$pos,
               as.integer(x@other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1))
  expect_equal(as.character(map$chrom),
               as.character(x@other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1))
})

test_that("explicit snp.pos/snp.chr take precedence over a pre-set slot", {
  # F1 fix applied: an explicitly supplied snp.pos/snp.chr overwrites a
  # populated @position/@chromosome slot instead of being silently ignored
  skip_if_not_installed("dartR.data")
  x <- vcf_fixture()
  x@position <- as.integer(seq(11, by = 11, length.out = nLoc(x)))
  x@chromosome <- factor(rep("77", nLoc(x)))
  map <- run_gl2vcf_map(x, "preset",
                        snp.pos = "ChromPos_Platypus_Chrom_NCBIv1",
                        snp.chr = "Chrom_Platypus_Chrom_NCBIv1")
  expect_equal(map$pos,
               as.integer(x@other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1))
  expect_equal(as.character(map$chrom),
               as.character(x@other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1))
})

test_that("a valid pre-set position slot is used when no snp.pos is supplied", {
  # F1 fix applied: with no explicit arguments, valid slots are used as found
  skip_if_not_installed("dartR.data")
  x <- vcf_fixture()
  x@position <- as.integer(seq(11, by = 11, length.out = nLoc(x)))
  x@chromosome <- factor(rep("77", nLoc(x)))
  map <- run_gl2vcf_map(x, "slotonly")
  expect_equal(map$pos, seq(11L, by = 11L, length.out = 8))
  expect_equal(as.character(map$chrom), rep("77", 8))
})

test_that("a factor-typed snp.pos field yields the labelled positions", {
  # F5 fix applied: coercion via as.integer(as.character(...)) returns the
  # factor labels, not the level codes
  skip_if_not_installed("dartR.data")
  x <- vcf_fixture()
  x@other$loc.metrics$posfac <-
    factor(x@other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1)
  map <- run_gl2vcf_map(x, "factor", snp.pos = "posfac",
                        snp.chr = "Chrom_Platypus_Chrom_NCBIv1")
  truth <- x@other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
  expect_equal(map$pos, as.integer(truth))                    # F5 fix applied
})

test_that("a non-numeric snp.pos field is fatal, not silently NA", {
  # F5 fix applied: values that cannot be interpreted as integer
  # positions stop with a clear message
  skip_if_not_installed("dartR.data")
  x <- vcf_fixture()
  x@other$loc.metrics$badpos <- rep("not-a-number", nLoc(x))
  td <- dummy_plink_dir("badpos")
  expect_error(
    suppressMessages(capture.output(
      gl2vcf(x, plink.bin.path = td, outfile = "o", outpath = td,
             snp.pos = "badpos", verbose = 0))),
    "cannot be interpreted as integer SNP positions")
})

test_that("a wrong-length position slot is zeroed with a gated warning", {
  skip_if_not_installed("dartR.data")
  x <- vcf_fixture()
  x@position <- c(1000L, 2000L)
  map <- run_gl2vcf_map(x, "wronglen")
  expect_equal(map$pos, rep(0L, 8))
  # F7 fix applied: discarding a malformed slot warns at verbose >= 1
  td <- dummy_plink_dir("wronglenv1")
  o <- capture.output(try(suppressMessages(
    gl2vcf(x, plink.bin.path = td, outfile = "out", outpath = td,
           verbose = 1)), silent = TRUE))
  expect_true(any(grepl("@position slot has length", o)))
})

test_that("a nonexistent snp.pos field is fatal with a clear message", {
  skip_if_not_installed("dartR.data")
  td <- dummy_plink_dir("nofield")
  expect_error(
    suppressMessages(capture.output(
      gl2vcf(vcf_fixture(), plink.bin.path = td, outfile = "o",
             outpath = td, snp.pos = "no_such_field", verbose = 0))),
    "not present in loc.metrics")
})

test_that("SilicoDArT is rejected explicitly", {
  # F4 fix applied: accept = 'SNP' (DAT7) makes presence/absence data a
  # clear datatype error instead of a cryptic failure deep in gl2plink
  td <- file.path(tempdir(), "vcfgs")
  dir.create(td, showWarnings = FALSE)
  expect_error(
    suppressMessages(capture.output(
      gl2vcf(testset2.gs[1:4, 1:6], plink.bin.path = td, outfile = "o",
             outpath = td, verbose = 0))),
    "SilicoDArT")                                             # F4 fix applied
})

test_that("a missing PLINK binary is fatal up front with the download URL", {
  # F6 fix applied: the binary is checked before any work is done
  skip_if_not_installed("dartR.data")
  td <- file.path(tempdir(), "vcfnoplink")
  dir.create(td, showWarnings = FALSE)
  expect_error(
    suppressMessages(capture.output(
      gl2vcf(vcf_fixture(), plink.bin.path = td, outfile = "o",
             outpath = td, verbose = 0))),
    "cog-genomics")
})

test_that("end-to-end VCF via PLINK: coordinates, alleles, genotypes", {
  skip_if_not_installed("dartR.data")
  plink_dir <- Sys.getenv("PLINK19_DIR")
  skip_if(plink_dir == "" ||
            !any(file.exists(file.path(plink_dir, c("plink", "plink.exe")))),
          "PLINK 1.9 not available (set PLINK19_DIR)")
  td <- file.path(tempdir(), "vcfe2e")
  dir.create(td, showWarnings = FALSE)
  x <- vcf_fixture()
  # F3 fix applied: verbose = 0 is silent -- gl2plink runs at the passed
  # verbosity and PLINK's captured log prints only at verbose >= 2
  o <- capture.output(m <- capture.output(
    gl2vcf(x, plink.bin.path = plink_dir, outfile = "e2e", outpath = td,
           snp.pos = "ChromPos_Platypus_Chrom_NCBIv1",
           snp.chr = "Chrom_Platypus_Chrom_NCBIv1", verbose = 0),
    type = "message"))
  expect_length(o, 0)                                         # F3 fix applied
  expect_length(m, 0)                                         # F3 fix applied
  v <- readLines(file.path(td, "e2e.vcf"))
  expect_equal(v[1], "##fileformat=VCFv4.2")
  body <- strsplit(v[!grepl("^#", v)], "\t")
  expect_length(body, nLoc(x))
  ids <- vapply(body, `[`, "", 3)
  expect_setequal(ids, locNames(x))
  # POS carries the genome coordinate for each locus
  pos <- as.integer(vapply(body, `[`, "", 2))
  truth <- x@other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
  expect_equal(pos, truth[match(ids, locNames(x))])
  # REF/ALT versus loc.all, and dosage -> GT for the first individual;
  # F2 fix applied: REF is pinned from loc.all via --a2-allele, so loci
  # whose reference allele is unobserved in the sample keep their
  # recorded reference base instead of degrading to N
  la <- strsplit(x@loc.all, "/")
  gm <- as.matrix(x)
  for (k in seq_along(body)) {
    j <- match(ids[k], locNames(x))
    ref <- body[[k]][4]; alt <- body[[k]][5]
    expect_equal(alt, la[[j]][2])
    expect_equal(ref, la[[j]][1])                             # F2 fix applied
    gt <- body[[k]][10]
    expected_gt <- c(`0` = "0/0", `1` = "0/1", `2` = "1/1")[
      as.character(gm[1, j])]
    if (is.na(gm[1, j])) expected_gt <- "./."
    expect_equal(gt, unname(expected_gt))
  }
  # F8 fix applied: intermediates live in tempdir() and, with the PLINK
  # by-products (.log, .nosex), are removed on success
  expect_false(any(c("gl_plink_temp.map", "gl_plink_temp.ped",
                     "e2e.log", "e2e.nosex") %in% list.files(td)))
})
