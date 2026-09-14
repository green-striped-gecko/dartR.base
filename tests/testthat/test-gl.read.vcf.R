# Characterization tests for gl.read.vcf. Captured pre-review at
# upstream/dev ddaed27, then updated for the approved Phase C fixes
# (findings F1-F6, approved 2026-09-05); each changed expectation is
# annotated with its finding id. Findings and approval record:
# function-review/reports/dartR.base/gl.read.vcf.md.

write_vcf <- function(dir, name, records,
                      info_hdr = character(0), samples = c("S1", "S2", "S3")) {
  f <- file.path(dir, name)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    "##contig=<ID=chr2,length=2000>",
    info_hdr,
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", samples), collapse = "\t"),
    records), f)
  f
}

vcf_tmpdir <- function() {
  d <- tempfile("vcffix")
  dir.create(d)
  d
}

info2 <- c(
  "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
  "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">")

basic_records <- c(
  "chr1\t100\trs1\tA\tG\t50\tPASS\tDP=10;AC=2\tGT\t0/0\t0/1\t1/1",
  "chr1\t200\trs2\tC\tT\t60\tPASS\tDP=12;AC=3\tGT\t0|1\t1|1\t0|0",
  "chr2\t150\trs3\tG\tA\t70\tPASS\tDP=8;AC=1\tGT\t./.\t0/0\t0/1")

test_that("gl.read.vcf reads biallelic SNPs: ALT dosage, slots, loc.metrics", {
  skip_if_not_installed("vcfR")
  d <- vcf_tmpdir()
  f <- write_vcf(d, "basic.vcf", basic_records, info_hdr = info2)

  capture.output(x <- suppressWarnings(gl.read.vcf(f, verbose = 0)))

  expect_s4_class(x, "dartR")
  expect_equal(c(nInd(x), nLoc(x)), c(3L, 3L))
  expect_equal(unique(as.integer(ploidy(x))), 2L)
  expect_equal(indNames(x), c("S1", "S2", "S3"))
  expect_equal(locNames(x), c("rs1", "rs2", "rs3"))

  # Dosage counts the ALT allele; phased separators handled; ./. becomes NA.
  expected <- matrix(
    c(0, 1, NA,
      1, 2, 0,
      2, 0, 1),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("S1", "S2", "S3"), c("rs1", "rs2", "rs3")))
  expect_equal(as.matrix(x), expected)

  # Genome coordinates land in the slots (PR #330 convention).
  expect_equal(x@position, c(100L, 200L, 150L))
  expect_equal(as.character(x@chromosome), c("chr1", "chr1", "chr2"))
  expect_equal(x@loc.all, c("A/G", "C/T", "G/A"))

  # QUAL/FILTER/INFO fields land in loc.metrics.
  lm <- x@other$loc.metrics
  expect_equal(as.numeric(lm$QUAL), c(50, 60, 70))
  expect_equal(as.character(lm$FILTER), c("PASS", "PASS", "PASS"))
  expect_equal(as.numeric(as.character(lm$DP)), c(10, 12, 8))
  expect_equal(as.numeric(as.character(lm$AC)), c(2, 3, 1))
  # history holds the gl.read.vcf call plus the entry appended by the
  # internal gl.recalc.metrics call.
  expect_equal(length(x@other$history), 2L)
})

test_that("gl.read.vcf is silent at verbose = 0", {
  skip_if_not_installed("vcfR")
  d <- vcf_tmpdir()
  f <- write_vcf(d, "basic.vcf", basic_records, info_hdr = info2)
  out <- capture.output(x <- suppressWarnings(gl.read.vcf(f, verbose = 0)))
  # [approved F3] Previously 48 lines printed: gl.compliance.check() and
  # gl.recalc.metrics() ran at their own default verbosity. The user's
  # verbose is now passed through.
  expect_length(out, 0)
})

test_that("gl.read.vcf drops multi-allelic records consistently and preserves FILTER", {
  skip_if_not_installed("vcfR")
  d <- vcf_tmpdir()
  recs <- c(
    "chr1\t100\trs1\tA\tG\t50\tPASS\tDP=10;AC=2\tGT\t0/0\t0/1\t1/1",
    "chr1\t200\trs2\tC\tG,T\t60\tPASS\tDP=12;AC=1,2\tGT\t1/2\t0/1\t0/2",
    "chr1\t300\trs3\tA\tAT\t70\tPASS\tDP=9;AC=2\tGT\t0/1\t1/1\t0/0",
    "chr1\t400\trs4\tG\tA\t80\tPASS\tDP=7;AC=1\tGT\t0/0\t0/1\t0/0")
  f <- write_vcf(d, "acmulti.vcf", recs, info_hdr = info2)

  capture.output(x <- suppressWarnings(gl.read.vcf(f, verbose = 0)))

  # rs2 (multi-allelic) dropped from genotypes AND from slots/metrics in step.
  expect_equal(locNames(x), c("rs1", "rs3", "rs4"))
  expect_equal(x@position, c(100L, 300L, 400L))
  expect_equal(length(x@chromosome), nLoc(x))
  expect_equal(nrow(x@other$loc.metrics), nLoc(x))
  # The indel rs3 (A/AT) is retained silently and coded as an ALT dosage.
  expect_equal(x@loc.all[2], "A/AT")
  expect_equal(unname(as.matrix(x)[, "rs3"]), c(1, 2, 0))
  # [approved F4] Previously every info column was coerced numeric in this
  # branch, so FILTER "PASS" became NA for every retained locus. Only
  # numeric-parseable columns are coerced now.
  expect_equal(as.character(x@other$loc.metrics$FILTER),
               c("PASS", "PASS", "PASS"))
  expect_equal(as.numeric(x@other$loc.metrics$DP), c(10, 9, 7))
})

test_that("gl.read.vcf assigns INFO values by key when key order varies across records", {
  skip_if_not_installed("vcfR")
  d <- vcf_tmpdir()
  recs <- c(
    "chr1\t100\tw1\tA\tG\t50\tPASS\tDP=10;AC=2\tGT\t0/0\t0/1\t1/1",
    "chr1\t200\tw2\tC\tT\t60\tPASS\tAC=3;DP=99\tGT\t0/1\t1/1\t0/0")
  f <- write_vcf(d, "infoswap.vcf", recs, info_hdr = info2)

  capture.output(x <- suppressWarnings(gl.read.vcf(f, verbose = 0)))

  # [approved F2] Previously column names came from record 1 and record
  # 2's values were filled by position, so w2 got DP=3 and AC=99. INFO is
  # now parsed per record by key: w2 has DP=99, AC=3 (the truth).
  lm <- x@other$loc.metrics
  expect_equal(as.numeric(as.character(lm$DP)), c(10, 99))
  expect_equal(as.numeric(as.character(lm$AC)), c(2, 3))
})

test_that("gl.read.vcf dosage mode keeps data-derived ploidy with polyploid dosages", {
  skip_if_not_installed("vcfR")
  d <- vcf_tmpdir()
  recs <- c(
    "chr1\t100\tt1\tA\tG\t.\tPASS\t.\tGT\t0/0/1\t1/1/1\t0/0/0",
    "chr1\t200\tt2\tC\tT\t.\tPASS\t.\tGT\t0/1/1\t0/0/0\t1/1/1",
    "chr1\t300\tt3\tG\tA\t.\tPASS\t.\tGT\t0/0/0\t0/0/1\t././.")
  f <- write_vcf(d, "poly.vcf", recs, samples = c("T1", "T2", "T3"))

  capture.output(x <- suppressWarnings(gl.read.vcf(f, mode = "dosage",
                                                   verbose = 0)))

  # Hand-computed ALT copy numbers, including 0/0/1 -> 1 and 1/1/1 -> 3.
  expected <- matrix(
    c(1, 2, 0,
      3, 0, 1,
      0, 3, NA),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("T1", "T2", "T3"), c("t1", "t2", "t3")))
  expect_equal(as.matrix(x), expected)
  # [approved F1] Previously every individual was stamped diploid
  # (ploidy 2,2,2) although the matrix carries dosages of 3. Dosage mode
  # now stamps the data's maximum copy number (uniform, as
  # gl.compliance.check requires a single ploidy level per object).
  expect_equal(as.integer(ploidy(x)), c(3L, 3L, 3L))
})

test_that("gl.read.vcf genotype mode recodes haploid calls as diploid homozygotes, warning at verbose >= 1", {
  skip_if_not_installed("vcfR")
  d <- vcf_tmpdir()
  recs <- c(
    "chr1\t100\th1\tA\tG\t.\tPASS\t.\tGT\t0\t1\t0",
    "chr1\t200\th2\tC\tT\t.\tPASS\t.\tGT\t1\t1\t0")
  f <- write_vcf(d, "haploid.vcf", recs, samples = c("H1", "H2", "H3"))

  capture.output(x <- suppressWarnings(gl.read.vcf(f, verbose = 0)))

  expected <- matrix(
    c(0, 2,
      2, 2,
      0, 0),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("H1", "H2", "H3"), c("h1", "h2")))
  expect_equal(as.matrix(x), expected)
  expect_equal(as.integer(ploidy(x)), c(2L, 2L, 2L))

  # [approved F6] The recoding is unchanged (documented diploid coding),
  # but it is no longer silent: a warning states it at verbose >= 1.
  out <- capture.output(x1 <- suppressWarnings(gl.read.vcf(f, verbose = 1)))
  expect_true(any(grepl("non-diploid genotype calls detected", out)))
})

test_that("gl.read.vcf ind.metafile realigns by id and retains unmatched individuals", {
  skip_if_not_installed("vcfR")
  d <- vcf_tmpdir()
  f <- write_vcf(d, "basic.vcf", basic_records, info_hdr = info2)
  im <- file.path(d, "indmeta.csv")
  # scrambled order, S2 deliberately absent
  write.csv(data.frame(id = c("S3", "S1"), pop = c("PZ", "PA"),
                       age = c(3, 1)), im, row.names = FALSE)

  capture.output(x <- suppressWarnings(
    gl.read.vcf(f, ind.metafile = im, verbose = 0)))

  # [approved F5] Previously S2 was silently dropped from the returned
  # object (2 individuals came back). Unmatched individuals are now
  # retained with NA metadata; matched rows are still realigned by id.
  expect_equal(indNames(x), c("S1", "S2", "S3"))
  expect_equal(as.character(x@other$ind.metrics$id), c("S1", "S2", "S3"))
  expect_equal(x@other$ind.metrics$age, c(1, NA, 3))
  expect_equal(as.character(pop(x))[c(1, 3)], c("PA", "PZ"))
  expect_true(is.na(as.character(pop(x))[2]))
  # Genotypes keep all three individuals in SNP-data order.
  expect_equal(unname(as.matrix(x)[, "rs1"]), c(0, 1, 2))

  # The drop-list warning prints at verbose >= 1.
  out <- capture.output(x1 <- suppressWarnings(
    gl.read.vcf(f, ind.metafile = im, verbose = 1)))
  expect_true(any(grepl("retained with NA metadata", out)))
})

test_that("gl.read.vcf reads a gzipped vcf", {
  skip_if_not_installed("vcfR")
  d <- vcf_tmpdir()
  f <- write_vcf(d, "basic.vcf", basic_records, info_hdr = info2)
  gz <- file.path(d, "basic.vcf.gz")
  con <- gzfile(gz, "wb"); writeLines(readLines(f), con); close(con)

  capture.output(x <- suppressWarnings(gl.read.vcf(gz, verbose = 0)))
  expect_equal(c(nInd(x), nLoc(x)), c(3L, 3L))
})

test_that("gl.read.vcf never prints its own FLAG SCRIPT END (current gap)", {
  skip_if_not_installed("vcfR")
  d <- vcf_tmpdir()
  f <- write_vcf(d, "basic.vcf", basic_records, info_hdr = info2)
  out <- capture.output(invisible(suppressWarnings(gl.read.vcf(f, verbose = 2))))
  expect_false(any(grepl("Completed: gl.read.vcf", out)))
})
