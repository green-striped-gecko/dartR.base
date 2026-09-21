test_that("gl.filter.contamination removes exactly the individuals the report flags", {
  x <- testset.gl
  capture.output(r <- gl.report.contamination(x, plot.display = FALSE, verbose = 0))
  expected.drop <- r$ind$id[r$ind$flag %in% c("suspect", "adjacent")]
  out <- capture.output(x2 <- gl.filter.contamination(x, verbose = 0))
  expect_length(out, 0)
  expect_s4_class(x2, "genlight")
  expect_equal(nInd(x2), nInd(x) - length(expected.drop))
  expect_setequal(indNames(x2), setdiff(indNames(x), expected.drop))
  expect_equal(nLoc(x2), nLoc(x))
  expect_true(all(ploidy(x2) == 2))
  expect_equal(nrow(x2@other$ind.metrics), nInd(x2))
  expect_length(x2@other$history, length(x@other$history) + 1)
})

test_that("a synthetic mixture of two individuals is removed", {
  # Same construction as the report test: the host is the WA animal with
  # median heterozygosity, the donor the first SA animal
  set.seed(1)
  x <- bandicoot.gl
  m <- as.matrix(x)
  het0 <- rowMeans(m == 1, na.rm = TRUE)
  wa <- which(pop(x) == "WA")
  host <- wa[which.min(abs(het0[wa] - median(het0[wa])))]
  donor <- which(pop(x) == "SA")[1]
  loci <- which(!is.na(m[host, ]) & !is.na(m[donor, ]) & m[host, ] != m[donor, ])
  m[host, sample(loci, round(0.5 * length(loci)))] <- 1
  x2 <- x
  x2@gen <- new("genlight", m, ploidy = 2)@gen
  indNames(x2) <- indNames(x); pop(x2) <- pop(x)
  out <- capture.output(x3 <- gl.filter.contamination(x2, verbose = 3))
  expect_false(indNames(x)[host] %in% indNames(x3))
  expect_true(indNames(x)[donor] %in% indNames(x3))
  expect_true(any(grepl(indNames(x)[host], out, fixed = TRUE)))
  expect_true(any(grepl("Number of individuals removed:", out)))
})

test_that("flag selects the classes removed", {
  x <- testset.gl
  capture.output(r <- gl.report.contamination(x, plot.display = FALSE, verbose = 0))
  # only plate-adjacent suspects: testset.gl has no plate positions, so none
  capture.output(x.adj <- gl.filter.contamination(x, flag = "adjacent", verbose = 0))
  expect_equal(nInd(x.adj), nInd(x))
  # all three classes
  capture.output(x.all <- gl.filter.contamination(
    x, flag = c("suspect", "adjacent", "rare-only"), verbose = 0))
  expect_equal(nInd(x.all), sum(r$ind$flag == ""))
  expect_error(capture.output(gl.filter.contamination(x, flag = "bogus", verbose = 0)),
               "flag must be one or more of")
})

test_that("locus metric flags are reset when individuals are removed, or recalculated", {
  x <- testset.gl
  capture.output(r <- gl.report.contamination(x, plot.display = FALSE, verbose = 0))
  skip_if(sum(r$ind$flag %in% c("suspect", "adjacent")) == 0,
          "testset.gl yields no suspect at the default thresholds")
  capture.output(x2 <- gl.filter.contamination(x, verbose = 0))
  expect_false(isTRUE(x2@other$loc.metrics.flags$CallRate))
  expect_false(isTRUE(x2@other$loc.metrics.flags$maf))
  capture.output(x3 <- gl.filter.contamination(x, recalc = TRUE, verbose = 0))
  expect_true(isTRUE(x3@other$loc.metrics.flags$CallRate))
  expect_equal(nrow(x3@other$loc.metrics), nLoc(x3))
  capture.output(x4 <- gl.filter.contamination(x, mono.rm = TRUE, verbose = 0))
  expect_true(nLoc(x4) <= nLoc(x))
  expect_equal(nrow(x4@other$loc.metrics), nLoc(x4))
})

test_that("presence/absence data and bad screen parameters stop", {
  expect_error(capture.output(gl.filter.contamination(testset.gs, verbose = 0)))
  expect_error(capture.output(gl.filter.contamination(testset.gl, z.flag = 0, verbose = 0)))
  expect_error(capture.output(gl.filter.contamination(testset.gl, plate = data.frame(id = 1),
                                                       verbose = 0)))
})
