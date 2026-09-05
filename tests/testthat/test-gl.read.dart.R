# Baseline characterization tests for gl.read.dart (Phase A, dartr-function-review).
# These snapshot CURRENT behaviour on the reviewed commit, including known
# defects, so that any change in behaviour is detected. They are not an
# endorsement of that behaviour; see
# function-review/reports/dartR.base/gl.read.dart.md.
#
# NOTE on @position: open PR #330 (position-genome-only) changes @position/
# @chromosome to NULL after import; pre-#330 code fills @position with a copy
# of loc.metrics$SnpPosition. This branch is based on pre-#330 dev, so the
# position expectations below accept either state coherently and stay valid
# when #330 lands.

dartfile <- system.file("extdata", "testset_SNPs_2Row.csv", package = "dartR.data")
metafile <- system.file("extdata", "testset_metadata.csv", package = "dartR.data")

# Read once, reuse across tests
out <- capture.output(suppressWarnings(
  gl <- gl.read.dart(dartfile, ind.metafile = metafile, verbose = 0)
))

test_that("gl.read.dart returns a compliant dartR SNP object from the 2-row testset", {
  expect_s4_class(gl, "dartR")
  expect_equal(nInd(gl), 250)
  expect_equal(nLoc(gl), 255)
  expect_true(all(ploidy(gl) == 2))
  expect_equal(nPop(gl), 30)
  expect_length(gl@other$history, 3)
})

test_that("@position tracks loc.metrics$SnpPosition (or is NULL post-#330)", {
  # [#330 interaction] pre-#330: @position holds SnpPosition; post-#330: NULL
  expect_true(is.null(gl@position) ||
                all(as.integer(gl@position) ==
                      as.integer(gl@other$loc.metrics$SnpPosition)))
  expect_true("SnpPosition" %in% names(gl@other$loc.metrics))
  expect_equal(head(gl@other$loc.metrics$SnpPosition, 5), c(19, 38, 24, 33, 21))
})

test_that("loc.metrics rows register 1:1 with loci and carry derived metrics", {
  expect_equal(nrow(gl@other$loc.metrics), nLoc(gl))
  expect_true(all(c("AlleleID", "CloneID", "SNP", "SnpPosition", "CallRate",
                    "RepAvg", "clone", "uid", "rdepth", "maf") %in%
                    names(gl@other$loc.metrics)))
  # maf is present in the final object; the utils.recalc.maf result is now
  # assigned inside gl.read.dart # [approved F2]
  expect_false(any(is.na(gl@other$loc.metrics$maf)))
  expect_true(gl@other$loc.metrics.flags$maf)
  # SilicoDArT-only flags are off for SNP data
  expect_false(gl@other$loc.metrics.flags$OneRatio)
  expect_false(gl@other$loc.metrics.flags$PIC)
})

test_that("2-row genotype encoding matches hand-read csv cells", {
  # From testset_SNPs_2Row.csv (topskip = 3, pairs ref-row/snp-row):
  # AA010915: locus1 0/1 -> 2 (hom alt), locus2 -/- -> NA,
  #           locus3 1/0 -> 0 (hom ref), locus4 1/0 -> 0
  # UC_00126: same four raw pairs
  m <- as.matrix(gl)
  expect_equal(unname(m["AA010915", 1:4]), c(2, NA, 0, 0))
  expect_equal(unname(m["UC_00126", 1:4]), c(2, NA, 0, 0))
})

test_that("ind.metrics carry metafile columns plus service/plate columns", {
  expect_true(all(c("id", "pop", "lat", "lon", "sex", "maturity",
                    "service", "plate_location") %in% names(gl@other$ind.metrics)))
  expect_equal(nrow(gl@other$ind.metrics), nInd(gl))
  # CHARACTERIZES A DEFECT: with only 3 header rows in this file, the default
  # service.row = 1 / plate.row = 3 read past the header block, so 'service'
  # holds the plate-well row and plate_location is built from the header and
  # first data row (e.g. "UC_1-AA0109150"). See report finding on bounds check.
  expect_equal(as.character(gl@other$ind.metrics$service[1]), "A5")
  expect_equal(as.character(gl@other$ind.metrics$plate_location[1]), "UC_1-AA0109150")
})

test_that("verbose = 0 is silent: gl.compliance.check log no longer leaks", {
  # gl.read.dart now passes verbose through to gl.compliance.check, so a
  # verbose = 0 read prints nothing # [approved F1]
  expect_length(out, 0)
})

test_that("explicit topskip skips service/plate extraction (values all NA)", {
  o <- capture.output(suppressWarnings(
    gl3 <- gl.read.dart(dartfile, ind.metafile = metafile, topskip = 3, verbose = 0)
  ))
  expect_equal(nInd(gl3), 250)
  expect_equal(nLoc(gl3), 255)
  expect_true(all(is.na(gl3@other$ind.metrics$service)))
  expect_true(all(is.na(gl3@other$ind.metrics$plate_location)))
})

test_that("1-row format: autodetection, encoding, and service/plate extraction", {
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
  o <- capture.output(suppressWarnings(g1 <- gl.read.dart(f, verbose = 0)))
  expect_equal(nInd(g1), 3)
  expect_equal(nLoc(g1), 3)
  expect_true(all(ploidy(g1) == 2))
  m <- as.matrix(g1)
  # 1-row coding: raw 0 = hom ref -> 0; raw 1 = hom SNP -> 2; raw 2 = het -> 1
  # raw columns per individual: ind1 (0,2,1); ind2 (1,1,NA); ind3 (2,0,0)
  expect_equal(unname(m["ind1", ]), c(0, 1, 2))
  expect_equal(unname(m["ind2", ]), c(2, 2, NA))
  expect_equal(unname(m["ind3", ]), c(1, 0, 0))
  # with 5 header rows the default service/plate rows are in range and correct
  expect_equal(as.character(g1@other$ind.metrics$service),
               c("ServiceA", "ServiceB", "ServiceC"))
  expect_equal(as.character(g1@other$ind.metrics$plate_location),
               c("P1-A1", "P1-B2", "P1-C3"))
  # [#330 interaction] pre-#330: @position holds SnpPosition; post-#330: NULL
  expect_true(is.null(g1@position) ||
                all(as.integer(g1@position) ==
                      as.integer(g1@other$loc.metrics$SnpPosition)))
})

test_that("file without '*' rows fails with a clear lastmetric error", {
  # lastmetric autodetection now runs after the preamble and reports the
  # missing '*' block instead of "argument of length 0" # [approved F3]
  f <- tempfile(fileext = ".csv")
  writeLines(c(
    "AlleleID,SNP,SnpPosition,CallRate,AvgCountRef,AvgCountSnp,RepAvg,ind1,ind2,ind3",
    "100001|F|0-10:A>G,10:A>G,10,1,20,15,1,0,1,2",
    "100002|F|0-5:C>T,5:C>T,5,1,18,12,0.98,2,1,0"), f)
  expect_error(suppressWarnings(gl.read.dart(f, topskip = 0, verbose = 0)),
               "lastmetric parameter")
})

test_that("a mistyped filename fails with a clear error", {
  # file existence is now checked before any read # [approved F3]
  expect_error(suppressWarnings(gl.read.dart("no_such_file.csv", verbose = 0)),
               "does not exist")
})

test_that("missing AvgCount columns leave rdepth NA instead of crashing", {
  # the rdepth block is now guarded on CallRate/AvgCountRef/AvgCountSnp being
  # present; when absent, rdepth is set NA with a gated warning
  # # [approved F4]
  f <- tempfile(fileext = ".csv")
  writeLines(c(
    "*,*,*,*,*,ServiceA,ServiceB,ServiceC",
    "*,*,*,*,*,BC1,BC2,BC3",
    "*,*,*,*,*,P1,P1,P1",
    "*,*,*,*,*,A,B,C",
    "*,*,*,*,*,1,2,3",
    "AlleleID,SNP,SnpPosition,CallRate,RepAvg,ind1,ind2,ind3",
    "100001|F|0-10:A>G,10:A>G,10,1,1,0,1,2",
    "100002|F|0-5:C>T,5:C>T,5,1,0.98,2,1,0"), f)
  o <- capture.output(suppressWarnings(g2 <- gl.read.dart(f, verbose = 0)))
  expect_equal(nInd(g2), 3)
  expect_equal(nLoc(g2), 2)
  expect_true(all(is.na(g2@other$loc.metrics$rdepth)))
})
