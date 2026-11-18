# Characterization tests for utils.vcfr2genlight.polyploid. Captured
# pre-review at upstream/dev ddaed27, then updated for the approved
# Phase C fixes (findings F1-F3, approved 2026-09-05); each changed
# expectation is annotated with its finding id. Findings and approval
# record: function-review/reports/dartR.base/utils.vcfr2genlight.polyploid.md.

read_vcfr <- function(dir, name, records, samples) {
  f <- file.path(dir, name)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", samples), collapse = "\t"),
    records), f)
  vcfR::read.vcfR(f, verbose = FALSE)
}

vcfr_tmpdir <- function() {
  d <- tempfile("vcfrfix")
  dir.create(d)
  d
}

poly_records <- c(
  "chr1\t100\tt1\tA\tG\t.\tPASS\t.\tGT\t0/0/1\t1/1/1\t0/0/0",
  "chr1\t200\tt2\tC\tT\t.\tPASS\t.\tGT\t0/1/1\t0/0/0\t1/1/1",
  "chr1\t300\tt3\tG\tA\t.\tPASS\t.\tGT\t0/0/0\t0/0/1\t././.")

test_that("polyploid dosage mode counts ALT copies; ploidy from GT arity", {
  skip_if_not_installed("vcfR")
  d <- vcfr_tmpdir()
  v <- read_vcfr(d, "poly.vcf", poly_records, c("T1", "T2", "T3"))

  x <- suppressWarnings(utils.vcfr2genlight.polyploid(x = v, mode2 = "dosage"))

  # Hand-computed ALT copy numbers for the triploid calls, 0/0/1 included.
  expected <- matrix(
    c(1, 2, 0,
      3, 0, 1,
      0, 3, NA),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("T1", "T2", "T3"), c("t1", "t2", "t3")))
  expect_equal(as.matrix(x), expected)

  # [approved F2] Ploidy is now the allele count (arity) of the GT calls:
  # 3,3,3 for uniformly triploid data. Previously adegenet inferred it
  # from each individual's max dosage, stamping T1 (no 1/1/1 call)
  # diploid: 2,3,3.
  expect_equal(as.integer(ploidy(x)), c(3L, 3L, 3L))

  # CHROM/POS/ID land in slots/locNames.
  expect_equal(as.character(x@position), c("100", "200", "300"))
  expect_equal(as.character(x@chromosome), c("chr1", "chr1", "chr1"))
  expect_equal(locNames(x), c("t1", "t2", "t3"))
})

test_that("polyploid genotype mode codes any mixed call as 1 and all-ALT as 2", {
  skip_if_not_installed("vcfR")
  d <- vcfr_tmpdir()
  v <- read_vcfr(d, "poly.vcf", poly_records, c("T1", "T2", "T3"))

  x <- suppressWarnings(utils.vcfr2genlight.polyploid(x = v, mode2 = "genotype"))

  expected <- matrix(
    c(1, 1, 0,
      2, 0, 1,
      0, 2, NA),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("T1", "T2", "T3"), c("t1", "t2", "t3")))
  expect_equal(as.matrix(x), expected)
  # [approved F2] Ploidy from GT arity: 3,3,3. Previously max-value
  # inference returned a haploid/diploid mixture (1,2,2) on this
  # uniformly triploid fixture.
  expect_equal(as.integer(ploidy(x)), c(3L, 3L, 3L))
})

test_that("half-missing calls are NA in both modes", {
  skip_if_not_installed("vcfR")
  d <- vcfr_tmpdir()
  recs <- c(
    "chr1\t100\tq1\tA\tG\t.\tPASS\t.\tGT\t./1\t0/0",
    "chr1\t200\tq2\tC\tT\t.\tPASS\t.\tGT\t0/.\t1/1")
  v <- read_vcfr(d, "partial.vcf", recs, c("P1", "P2"))

  # [approved F1] Previously genotype mode coded 0/. (no ALT allele
  # observed) and ./1 as heterozygous (1), and dosage mode returned 1 and
  # 0 for the same calls. Any call still carrying "." after separator
  # stripping is now missing data in both modes.
  x <- suppressWarnings(utils.vcfr2genlight.polyploid(x = v, mode2 = "genotype"))
  m <- as.matrix(x)
  expect_true(is.na(m["P1", "q1"]))
  expect_true(is.na(m["P1", "q2"]))
  expect_equal(unname(m["P2", ]), c(0, 2))

  xd <- suppressWarnings(utils.vcfr2genlight.polyploid(x = v, mode2 = "dosage"))
  md <- as.matrix(xd)
  expect_true(is.na(md["P1", "q1"]))
  expect_true(is.na(md["P1", "q2"]))
  expect_equal(unname(md["P2", ]), c(0, 2))

  # Ploidy from arity on the diploid calls.
  expect_equal(as.integer(ploidy(x)), c(2L, 2L))
})

test_that("multi-allelic loci are dropped with a warning before conversion", {
  skip_if_not_installed("vcfR")
  d <- vcfr_tmpdir()
  recs <- c(
    "chr1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1",
    "chr1\t200\trs2\tC\tG,T\t.\tPASS\t.\tGT\t1/2\t0/1")
  v <- read_vcfr(d, "multi.vcf", recs, c("S1", "S2"))

  expect_warning(x <- utils.vcfr2genlight.polyploid(x = v, mode2 = "genotype"),
                 "more than two alleles")
  expect_equal(locNames(x), "rs1")
  expect_equal(as.character(x@position), "100")
})

test_that("omitting mode2 defaults to genotype mode", {
  skip_if_not_installed("vcfR")
  d <- vcfr_tmpdir()
  recs <- "chr1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1"
  v <- read_vcfr(d, "one.vcf", recs, c("S1", "S2"))

  # [approved F3] Previously the signature default mode2 = mode resolved
  # to base::mode (a function) and the call failed with "comparison (==)
  # is possible only for atomic and list types". The default is now
  # "genotype", matching the caller.
  x <- utils.vcfr2genlight.polyploid(x = v)
  xg <- utils.vcfr2genlight.polyploid(x = v, mode2 = "genotype")
  expect_equal(as.matrix(x), as.matrix(xg))
  expect_equal(unname(as.matrix(x)[, 1]), c(0, 1))
})

test_that("an invalid mode2 stops with an empty condition message (current behaviour)", {
  skip_if_not_installed("vcfR")
  d <- vcfr_tmpdir()
  recs <- "chr1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1"
  v <- read_vcfr(d, "one.vcf", recs, c("S1", "S2"))

  # The message is cat()ed; stop() itself carries no message. (F4, LOW --
  # deferred, snapshot unchanged.)
  err <- tryCatch(
    capture.output(utils.vcfr2genlight.polyploid(x = v, mode2 = "bogus")),
    error = function(e) e)
  expect_s3_class(err, "error")
  expect_equal(conditionMessage(err), "")
})
