# Baseline characterization tests for utils.dart2genlight (Phase A,
# dartr-function-review). These snapshot CURRENT behaviour, including known
# defects; see function-review/reports/dartR.base/utils.dart2genlight.md.
#
# NOTE on @position: open PR #330 (position-genome-only) stops this function
# copying SnpPosition into @position (NULL thereafter); pre-#330 code fills
# @position with SnpPosition. This branch is based on pre-#330 dev, so the
# position expectation below accepts either state coherently and stays valid
# when #330 lands.

dartfile <- system.file("extdata", "testset_SNPs_2Row.csv", package = "dartR.data")
metafile <- system.file("extdata", "testset_metadata.csv", package = "dartR.data")

o <- capture.output(suppressWarnings(
  dart <- utils.read.dart(dartfile, verbose = 0)
))

test_that("2-row conversion: object shape, ploidy, encoding, metadata registration", {
  o2 <- capture.output(suppressWarnings(
    gl <- utils.dart2genlight(dart, ind.metafile = metafile,
                              probar = FALSE, verbose = 0)
  ))
  expect_s4_class(gl, "dartR")
  expect_equal(nInd(gl), 250)
  expect_equal(nLoc(gl), 255)
  expect_true(all(ploidy(gl) == 2))
  # [#330 interaction] pre-#330: @position holds SnpPosition; post-#330: NULL
  expect_true(is.null(gl@position) ||
                all(as.integer(gl@position) ==
                      as.integer(gl@other$loc.metrics$SnpPosition)))
  expect_null(gl@chromosome)
  expect_equal(nrow(gl@other$loc.metrics), nLoc(gl))
  expect_equal(nrow(gl@other$ind.metrics), nInd(gl))
  expect_equal(nPop(gl), 30)
  # hand-read cells (see test-gl.read.dart.R for the raw pairs)
  m <- as.matrix(gl)
  expect_equal(unname(m["AA010915", 1:4]), c(2, NA, 0, 0))
  expect_equal(unname(m["UC_00126", 1:4]), c(2, NA, 0, 0))
  # locus names are uid + allele pair from the SNP column
  expect_equal(locNames(gl)[1], "100049687-12-C/T")
})

test_that("individuals absent from ind.metafile are dropped, loudly at verbose >= 1", {
  # The drop is the approved contract (F1: explicit loud drop rather than
  # join semantics): dimensions are unchanged from baseline, but the removal
  # is now documented and announced with a REMOVED warning at verbose >= 1,
  # and verbose = 0 is silent # [approved F1, F5]
  sub <- read.csv(metafile)[1:100, ]
  f <- tempfile(fileext = ".csv")
  write.csv(sub, f, row.names = FALSE)
  o2 <- capture.output(suppressWarnings(
    gl <- utils.dart2genlight(dart, ind.metafile = f, probar = FALSE, verbose = 0)
  ))
  expect_equal(nInd(gl), 100)
  expect_equal(nrow(gl@other$ind.metrics), 100)
  expect_length(o2, 0)  # [approved F5]
  o3 <- capture.output(suppressWarnings(
    gl2 <- utils.dart2genlight(dart, ind.metafile = f, probar = FALSE, verbose = 1)
  ))
  expect_true(any(grepl("REMOVED", o3)))  # [approved F1]
})

test_that("metafile without a pop column assigns pop1 to everyone", {
  nopop <- read.csv(metafile)
  nopop$pop <- NULL
  f <- tempfile(fileext = ".csv")
  write.csv(nopop, f, row.names = FALSE)
  o2 <- capture.output(suppressWarnings(
    gl <- utils.dart2genlight(dart, ind.metafile = f, probar = FALSE, verbose = 0)
  ))
  expect_equal(nPop(gl), 1)
  expect_equal(popNames(gl), "pop1")
})

test_that("no metafile: single pop1 population, service/plate ind.metrics only", {
  o2 <- capture.output(suppressWarnings(
    gl <- utils.dart2genlight(dart, probar = FALSE, verbose = 0)
  ))
  expect_equal(nInd(gl), 250)
  expect_equal(popNames(gl), "pop1")
  expect_true(all(c("service", "plate_location") %in% names(gl@other$ind.metrics)))
})

test_that("a dart list without SNP/Variant columns fails with a clear header error", {
  # The fail-fast header check is reinstated: a missing SNP/Variant column
  # now raises a clear fatal error naming the expected columns instead of
  # "cannot coerce type 'closure'" # [approved F2]
  dart2 <- dart
  names(dart2$covmetrics)[names(dart2$covmetrics) == "SNP"] <- "SNPvariant"
  expect_error(
    suppressWarnings(capture.output(
      utils.dart2genlight(dart2, probar = FALSE, verbose = 0))),
    "'SNP' or 'Variant' column")
})

test_that("a 0/0 two-row pair converts to NA deliberately, without coercion warnings", {
  # '0/0' (double null) and partial-NA pairs are now explicit NA cases in
  # the translation table, so no "NAs introduced by coercion" warning is
  # emitted; genotype values are unchanged from baseline # [approved F3]
  f <- tempfile(fileext = ".csv")
  writeLines(c(
    "*,*,*,*,*,*,*,S1,S2,S3",
    "*,*,*,*,*,*,*,B1,B2,B3",
    "*,*,*,*,*,*,*,P1,P1,P1",
    "*,*,*,*,*,*,*,A,B,C",
    "*,*,*,*,*,*,*,1,2,3",
    "AlleleID,SNP,SnpPosition,CallRate,AvgCountRef,AvgCountSnp,RepAvg,ind1,ind2,ind3",
    "100001|F|0--10:A>G,,10,1,20,15,1,1,1,0",
    "100001|F|0-10:A>G-10:A>G,10:A>G,10,1,20,15,1,0,1,1",
    "100002|F|0--5:C>T,,5,1,18,12,0.98,0,0,1",
    "100002|F|0-5:C>T-5:C>T,5:C>T,5,1,18,12,0.98,0,1,0"), f)
  o2 <- capture.output(suppressWarnings(d <- utils.read.dart(f, verbose = 0)))
  expect_no_warning(
    capture.output(gl <- utils.dart2genlight(d, probar = FALSE, verbose = 0))
  )  # [approved F3]
  m <- as.matrix(gl)
  expect_true(is.na(m["ind1", "100002-5-C/T"]))
  expect_equal(unname(m["ind2", ]), c(1, 2))  # 1/1 -> 1 (het); 0/1 -> 2
})
