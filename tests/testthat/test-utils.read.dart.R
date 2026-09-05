# Baseline characterization tests for utils.read.dart (Phase A,
# dartr-function-review). These snapshot CURRENT behaviour, including known
# defects; see function-review/reports/dartR.base/utils.read.dart.md.

make_2row <- function(rows) {
  f <- tempfile(fileext = ".csv")
  writeLines(c(
    "*,*,*,*,*,*,*,ServiceA,ServiceB,ServiceC",
    "*,*,*,*,*,*,*,BC1,BC2,BC3",
    "*,*,*,*,*,*,*,P1,P1,P1",
    "*,*,*,*,*,*,*,A,B,C",
    "*,*,*,*,*,*,*,1,2,3",
    "AlleleID,SNP,SnpPosition,CallRate,AvgCountRef,AvgCountSnp,RepAvg,ind1,ind2,ind3",
    rows), f)
  f
}

test_that("real 2-row testset: dims, format, metrics, and header parsing", {
  dartfile <- system.file("extdata", "testset_SNPs_2Row.csv", package = "dartR.data")
  o <- capture.output(suppressWarnings(d <- utils.read.dart(dartfile, verbose = 0)))
  expect_equal(d$nrows, 2)
  expect_equal(d$nind, 250)
  expect_equal(d$nsnp, 255)
  expect_equal(nrow(d$covmetrics), 510)   # two rows per locus
  expect_equal(ncol(d$gendata), 250)
  expect_true(all(c("clone", "uid") %in% names(d$covmetrics)))
  # This file has only 3 header rows. service.row = 1 is in range, so service
  # still reads header row 1 (which in this file carries the plate wells,
  # e.g. "A5" - an in-range read that cannot be told apart from real service
  # codes). plate.row = 3 needs rows 3-5, which run past the header block, so
  # plate_location is now NA instead of a paste of the column-header and
  # first data rows # [approved F3]
  expect_equal(as.character(unlist(d$service))[1], "A5")
  expect_true(is.na(d$plate_location[1]))
})

test_that("true duplicate (quad) uid is removed and the file still parses", {
  f <- make_2row(c(
    "100001|F|0--10:A>G,,10,1,20,15,1,1,1,0",
    "100001|F|0-10:A>G-10:A>G,10:A>G,10,1,20,15,1,0,1,1",
    "100001|F|0--10:A>G,,10,1,20,15,1,1,1,0",
    "100001|F|0-10:A>G-10:A>G,10:A>G,10,1,20,15,1,0,1,1",
    "100002|F|0--5:C>T,,5,1,18,12,0.98,1,0,1",
    "100002|F|0-5:C>T-5:C>T,5:C>T,5,1,18,12,0.98,1,1,0",
    "100003|F|0--33:G>A,,33,1,25,10,1,1,1,1",
    "100003|F|0-33:G>A-33:G>A,33:G>A,33,1,25,10,1,0,0,1"))
  o <- capture.output(suppressWarnings(d <- utils.read.dart(f, verbose = 0)))
  expect_equal(d$nrows, 2)
  expect_equal(d$nsnp, 2)
  expect_equal(sort(unique(d$covmetrics$uid)), c("100002-5", "100003-33"))
})

test_that("an odd-count (triplet) uid is removed; the valid loci survive", {
  # Previously the duplicate-removal branch removed every uid with count > 1
  # (every valid 2-row locus) and the read died with a misleading format
  # error. Now only the uid whose row count departs from the modal count is
  # removed and the remaining loci parse # [approved F1]
  f <- make_2row(c(
    "100001|F|0--10:A>G,,10,1,20,15,1,1,1,0",
    "100001|F|0-10:A>G-10:A>G,10:A>G,10,1,20,15,1,0,1,1",
    "100001|F|0-10:A>G-10:A>G,10:A>G,10,1,20,15,1,0,1,1",
    "100002|F|0--5:C>T,,5,1,18,12,0.98,1,0,1",
    "100002|F|0-5:C>T-5:C>T,5:C>T,5,1,18,12,0.98,1,1,0",
    "100003|F|0--33:G>A,,33,1,25,10,1,1,1,1",
    "100003|F|0-33:G>A-33:G>A,33:G>A,33,1,25,10,1,0,0,1"))
  o <- capture.output(suppressWarnings(d <- utils.read.dart(f, verbose = 0)))
  expect_equal(d$nrows, 2)
  expect_equal(d$nsnp, 2)
  expect_equal(sort(unique(d$covmetrics$uid)), c("100002-5", "100003-33"))
})

test_that("1-row file with SNP codes 0/1/2 is detected as 1-row", {
  f <- tempfile(fileext = ".csv")
  writeLines(c(
    "*,*,*,*,*,*,*,ServiceA,ServiceB,ServiceC",
    "*,*,*,*,*,*,*,BC1,BC2,BC3",
    "*,*,*,*,*,*,*,P1,P1,P1",
    "*,*,*,*,*,*,*,A,B,C",
    "*,*,*,*,*,*,*,1,2,3",
    "AlleleID,SNP,SnpPosition,CallRate,AvgCountRef,AvgCountSnp,RepAvg,ind1,ind2,ind3",
    "100001|F|0-10:A>G,10:A>G,10,1,20,15,1,0,1,2",
    "100002|F|0-5:C>T,5:C>T,5,1,18,12,0.98,2,1,0",
    "100003|F|0-33:G>A,33:G>A,33,0.67,25,10,1,1,-,0"), f)
  o <- capture.output(suppressWarnings(d <- utils.read.dart(f, verbose = 0)))
  expect_equal(d$nrows, 1)
  expect_equal(d$nsnp, 3)
  expect_equal(as.character(unlist(d$service)), c("ServiceA", "ServiceB", "ServiceC"))
  expect_equal(d$plate_location, c("P1-A1", "P1-B2", "P1-C3"))
})

test_that("1-row file with no heterozygote codes now stops with a clear error", {
  # Previously format detection used max(genotype) only, so a 1-row file
  # scored 0/1 only (no hets called) was silently read as 2-row, halving the
  # locus count and pairing unrelated loci; only a console warning marked
  # it. The genotype-range and uid-count signals are now cross-checked and
  # the read stops when they disagree # [approved F2]
  f <- make_2row(c(
    "100001|F|0-10:A>G,10:A>G,10,1,20,15,1,0,1,0",
    "100002|F|0-5:C>T,5:C>T,5,1,18,12,0.98,1,1,0",
    "100003|F|0-33:G>A,33:G>A,33,1,25,10,1,1,0,0",
    "100004|F|0-7:T>C,7:T>C,7,1,22,11,1,0,0,1"))
  expect_error(
    suppressWarnings(capture.output(utils.read.dart(f, verbose = 0))),
    "contradicts that format")
})

test_that("duplicate individual names get the promised '_n' suffix", {
  # The make.unique(sep = '_') result is now applied to the gendata columns,
  # so the '_n' suffix promised by the warning is what downstream code sees
  # (previously read.csv's own '.1' suffix leaked through) # [approved F4]
  f <- tempfile(fileext = ".csv")
  writeLines(c(
    "*,*,*,*,*,*,*,ServiceA,ServiceB,ServiceC",
    "*,*,*,*,*,*,*,BC1,BC2,BC3",
    "*,*,*,*,*,*,*,P1,P1,P1",
    "*,*,*,*,*,*,*,A,B,C",
    "*,*,*,*,*,*,*,1,2,3",
    "AlleleID,SNP,SnpPosition,CallRate,AvgCountRef,AvgCountSnp,RepAvg,indA,indA,indB",
    "100001|F|0--10:A>G,,10,1,20,15,1,1,1,0",
    "100001|F|0-10:A>G-10:A>G,10:A>G,10,1,20,15,1,0,1,1",
    "100002|F|0--5:C>T,,5,1,18,12,0.98,0,0,1",
    "100002|F|0-5:C>T-5:C>T,5:C>T,5,1,18,12,0.98,0,1,0"), f)
  o <- capture.output(suppressWarnings(d <- utils.read.dart(f, verbose = 0)))
  expect_equal(colnames(d$gendata), c("indA", "indA_1", "indB"))  # [approved F4]
})

test_that("verbose = 0 is silent: all warnings are gated", {
  # The duplicate-individual warning is now gated at verbose >= 1, and no
  # other output survives verbose = 0 # [approved F5]
  f <- tempfile(fileext = ".csv")
  writeLines(c(
    "*,*,*,*,*,*,*,S1,S2",
    "*,*,*,*,*,*,*,B1,B2",
    "*,*,*,*,*,*,*,P1,P1",
    "*,*,*,*,*,*,*,A,B",
    "*,*,*,*,*,*,*,1,2",
    "AlleleID,SNP,SnpPosition,CallRate,AvgCountRef,AvgCountSnp,RepAvg,indA,indA",
    "100001|F|0--10:A>G,,10,1,20,15,1,1,1",
    "100001|F|0-10:A>G-10:A>G,10:A>G,10,1,20,15,1,0,1"), f)
  o <- capture.output(suppressWarnings(d <- utils.read.dart(f, verbose = 0)))
  expect_length(o, 0)  # [approved F5]
})
