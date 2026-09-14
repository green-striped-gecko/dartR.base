# Characterization tests for gl.allele.freq
# Baseline snapshotted before review (population-distance chain review,
# dev at ddaed27). Pins current behaviour, defects included; assertions
# tagged [pins defect] are expected to flip when the matching approved
# finding is applied.

test_that("popxloc cells match hand computation from as.matrix (SNP)", {
  glp <- gl.keep.pop(testset.gl, popNames(testset.gl)[1:4], verbose = 0)
  f <- gl.allele.freq(glp, percent = TRUE, by = "popxloc", verbose = 0)
  expect_s3_class(f, "data.frame")  # @return now documents a data.frame [approved F5]
  expect_equal(nrow(f), nPop(glp) * nLoc(glp))
  m <- as.matrix(glp)
  # hand: per pop x locus mean(genotype)/2 * 100, na.rm = TRUE, 2 dp
  for (p in popNames(glp)) {
    sub <- m[as.character(pop(glp)) == p, , drop = FALSE]
    hf <- round(colMeans(sub, na.rm = TRUE) / 2 * 100, 2)
    fr <- f[f$popn == p, ]
    fr <- fr[match(locNames(glp), as.character(fr$locus)), ]
    expect_equal(fr$frequency, unname(hf))
    expect_equal(fr$nobs, unname(colSums(!is.na(sub))))
    expect_equal(fr$nmissing, unname(colSums(is.na(sub))))
  }
})

test_that("locus 100049698-16-G/A anchors, including all-NA NaN cells", {
  glp <- gl.keep.pop(testset.gl, popNames(testset.gl)[1:4], verbose = 0)
  f <- gl.allele.freq(glp, percent = TRUE, by = "popxloc", verbose = 0)
  fr <- f[f$locus == "100049698-16-G/A", ]
  fr <- fr[order(as.character(fr$popn)), ]
  expect_equal(as.character(fr$popn),
               c("EmmacBrisWive", "EmmacBurdMist", "EmmacBurnBara",
                 "EmmacClarJack"))
  expect_equal(fr$frequency, c(0, NaN, 0, NaN))  # all-NA cells -> NaN
  expect_equal(fr$nmissing, c(3, 10, 8, 5))
})

test_that("percent=FALSE is the 2dp-rounded percentage divided by 100", {
  glp <- gl.keep.pop(testset.gl, popNames(testset.gl)[1:4], verbose = 0)
  fT <- gl.allele.freq(glp, percent = TRUE, by = "popxloc", verbose = 0)
  fF <- gl.allele.freq(glp, percent = FALSE, by = "popxloc", verbose = 0)
  expect_equal(fF$frequency, round(fT$frequency, 2) / 100)
})

test_that("by='pop' averages the popxloc cells across loci", {
  glp <- gl.keep.pop(testset.gl, popNames(testset.gl)[1:4], verbose = 0)
  fp <- gl.allele.freq(glp, percent = TRUE, by = "pop", verbose = 0)
  expect_equal(as.character(fp$popn),
               c("EmmacBrisWive", "EmmacBurdMist", "EmmacBurnBara",
                 "EmmacClarJack"))
  expect_equal(fp$frequency, c(35.7188, 34.8487, 35.5719, 34.3043))
  expect_equal(fp$nobs, c(8.7, 8.7, 9.6, 4.2))
})

test_that("by='loc' honours percent=TRUE", {
  # [approved F1] the by='loc' frequency column is now rescaled to the
  # percentage scale when percent = TRUE (previously always a proportion)
  fl <- gl.allele.freq(testset.gl, percent = TRUE, by = "loc", verbose = 0)
  expect_gt(max(fl$frequency, na.rm = TRUE), 1)
  expect_equal(fl$frequency,
               unname(round(colMeans(as.matrix(testset.gl),
                                     na.rm = TRUE) / 2 * 100, 4)))
  flF <- gl.allele.freq(testset.gl, percent = FALSE, by = "loc", verbose = 0)
  expect_equal(flF$frequency,
               unname(round(colMeans(as.matrix(testset.gl), na.rm = TRUE) / 2,
                            4)))
})

test_that("SilicoDArT: popxloc is presence percent but by='loc' is HALF", {
  gsp <- gl.keep.pop(testset.gs, popNames(testset.gs)[1:3], verbose = 0)
  m <- as.matrix(gsp)
  f <- gl.allele.freq(gsp, percent = TRUE, by = "popxloc", verbose = 0)
  p1 <- popNames(gsp)[1]
  sub <- m[as.character(pop(gsp)) == p1, , drop = FALSE]
  fr <- f[f$popn == p1, ]
  fr <- fr[match(locNames(gsp), as.character(fr$locus)), ]
  expect_equal(fr$frequency,
               unname(round(colMeans(sub, na.rm = TRUE) * 100, 2)))
  # [approved F2] by='loc' now returns the presence frequency itself
  # (previously the raw 0/1 matrix was divided by the SNP ploidy divisor,
  # halving it relative to the popxloc breakdown).
  fl <- gl.allele.freq(gsp, percent = FALSE, by = "loc", verbose = 0)
  expect_equal(fl$frequency,
               unname(round(colMeans(m, na.rm = TRUE), 4)))
})

test_that("simple=TRUE returns alf1/alf2 and silently overrides by/percent", {
  fs <- gl.allele.freq(testset.gl, simple = TRUE, verbose = 0)
  expect_equal(colnames(fs), c("alf1", "alf2"))
  expect_equal(nrow(fs), nLoc(testset.gl))
  expect_equal(fs$alf2,
               unname(round(colMeans(as.matrix(testset.gl),
                                     na.rm = TRUE) / 2, 4)))
  fso <- gl.allele.freq(testset.gl, simple = TRUE, percent = TRUE,
                        by = "pop", verbose = 0)
  expect_identical(fso, fs)  # override retained, now documented [approved F5]
})

test_that("an unrecognised 'by' stops with an informative error", {
  # [approved F4] by is validated against c('pop','loc','popxloc')
  expect_error(gl.allele.freq(testset.gl, by = "typo", verbose = 0),
               "popxloc")
})

test_that("verbose=0 is silent", {
  o <- capture.output(invisible(gl.allele.freq(testset.gl, verbose = 0)))
  expect_equal(length(o), 0)
})

test_that("a plain genlight without loc.metrics.flags is accepted", {
  # [approved F3] the monomorphs-flag access is now isTRUE()-guarded
  # (DAT5); a genlight not built by dartR runs instead of crashing with
  # "argument is of length zero"
  gg <- new("genlight", gen = as.matrix(testset.gl)[1:10, 1:20], ploidy = 2)
  pop(gg) <- pop(testset.gl)[1:10]
  f <- gl.allele.freq(gg, by = "popxloc", verbose = 0)
  expect_s3_class(f, "data.frame")
  expect_equal(nrow(f), nlevels(pop(gg)) * nLoc(gg))
})
