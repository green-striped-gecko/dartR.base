# Baseline characterization tests for gl.read.silicodart (Phase A,
# dartr-function-review). These snapshot CURRENT behaviour, including known
# defects; see function-review/reports/dartR.base/gl.read.silicodart.md.

silicofile <- system.file("extdata", "testset_SilicoDArT.csv", package = "dartR.data")
silicometa <- system.file("extdata", "testset_metadata_silicodart.csv", package = "dartR.data")

out <- capture.output(suppressWarnings(
  gs <- gl.read.silicodart(silicofile, ind.metafile = silicometa,
                           probar = FALSE, verbose = 0)
))

test_that("gl.read.silicodart returns a ploidy-1 dartR object from the testset", {
  expect_s4_class(gs, "dartR")
  expect_true(all(ploidy(gs) == 1))
  # the csv and the metafile both list the same 218 individuals, so no
  # drop occurs on the canonical fixture; the drop contract (explicit loud
  # drop, aligned with the SNP path) is exercised with a partial metafile
  # below # [approved F4]
  expect_equal(nInd(gs), 218)
  expect_equal(nLoc(gs), 255)
  expect_equal(nPop(gs), 29)
  m <- as.matrix(gs)
  expect_true(all(m %in% c(0, 1, NA)))
  expect_null(gs@position)
  expect_null(gs@chromosome)
})

test_that("loc.metrics register 1:1 with loci and keep Reproducibility", {
  expect_equal(nrow(gs@other$loc.metrics), nLoc(gs))
  expect_true(all(c("CloneID", "CallRate", "OneRatio", "PIC",
                    "Reproducibility") %in% names(gs@other$loc.metrics)))
  expect_true(all(c("id", "pop", "lat", "lon") %in% names(gs@other$ind.metrics)))
  expect_equal(nrow(gs@other$ind.metrics), nInd(gs))
})

test_that("verbose = 0 is silent: all messages gated, compliance check quiet", {
  # All messages are now verbosity-gated and gl.compliance.check receives
  # verbose, so a verbose = 0 read prints nothing # [approved F3]
  expect_length(out, 0)
})

test_that("metafile without a pop column assigns pop1 to everyone", {
  # The no-pop branch previously assigned pop(out) - an object that never
  # existed - and always crashed with "object 'out' not found". It now
  # mirrors the SNP path: every individual defaults to 'pop1'
  # # [approved F1]
  nopop <- read.csv(silicometa)
  nopop$pop <- NULL
  f <- tempfile(fileext = ".csv")
  write.csv(nopop, f, row.names = FALSE)
  o <- capture.output(suppressWarnings(
    gs2 <- gl.read.silicodart(silicofile, ind.metafile = f,
                              probar = FALSE, verbose = 0)
  ))
  expect_equal(nInd(gs2), 218)
  expect_equal(popNames(gs2), "pop1")
  expect_true(all(ploidy(gs2) == 1))
})

test_that("duplicate non-numeric CloneIDs get proper _n suffixed locus names", {
  # CloneID is converted to character before the uniquification loop, so the
  # pasted suffix strings are stored instead of factor-mismatch NAs
  # # [approved F2]
  f <- tempfile(fileext = ".csv")
  writeLines(c(
    "*,*,*,S1,S2",
    "*,*,*,B1,B2",
    "CloneID,AlleleSequence,Reproducibility,ind1,ind2",
    "CLONEA,ACGTACGT,1,1,0",
    "CLONEA,ACGTTCGT,0.95,0,1",
    "CLONEB,TTAAGGCC,1,1,1",
    "CLONEC,GGCCTTAA,0.99,0,1"), f)
  o <- capture.output(suppressWarnings(
    gL <- gl.read.silicodart(f, probar = FALSE, verbose = 0)
  ))
  expect_equal(sum(is.na(gL@other$loc.metrics$CloneID)), 0)  # [approved F2]
  expect_equal(sum(is.na(locNames(gL))), 0)                  # [approved F2]
  expect_true(all(c("CLONEA_1", "CLONEA_2") %in% locNames(gL)))
})

test_that("the partial-metafile drop is announced with a REMOVED warning at verbose 1", {
  # [approved F4] the drop itself is unchanged behaviour; it is now loud
  sub <- read.csv(silicometa)[1:100, ]
  f <- tempfile(fileext = ".csv")
  write.csv(sub, f, row.names = FALSE)
  o <- capture.output(suppressWarnings(
    gs3 <- gl.read.silicodart(silicofile, ind.metafile = f,
                              probar = FALSE, verbose = 1)
  ))
  expect_equal(nInd(gs3), 100)
  expect_true(any(grepl("REMOVED", o)))
})

test_that("an all-present matrix is rejected as fatal", {
  # CHARACTERIZES CURRENT BEHAVIOUR: min != 0 triggers the 0/1 coding check,
  # so a file where every tag is present in every individual cannot be read.
  f <- tempfile(fileext = ".csv")
  writeLines(c(
    "*,*,*,S1,S2",
    "*,*,*,B1,B2",
    "CloneID,AlleleSequence,Reproducibility,ind1,ind2",
    "CLONEA,ACGTACGT,1,1,1",
    "CLONEB,TTAAGGCC,1,1,1"), f)
  expect_error(
    suppressWarnings(capture.output(
      gl.read.silicodart(f, probar = FALSE, verbose = 0))),
    "must be 0 or 1")
})
