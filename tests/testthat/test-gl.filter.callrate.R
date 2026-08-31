test_that("gl.filter.callrate method='loc' removes loci below threshold (SNP)", {
  x <- testset.gl
  out <- capture.output(x2 <- gl.filter.callrate(x, method = "loc", threshold = 0.9,
                                                 plot.display = FALSE, verbose = 0))
  x_rc <- utils.recalc.callrate(x, verbose = 0)
  expected_keep <- sum(x_rc@other$loc.metrics$CallRate >= 0.9)
  expect_equal(nLoc(x2), expected_keep)
  expect_true(all(x2@other$loc.metrics$CallRate >= 0.9))
  expect_equal(nrow(x2@other$loc.metrics), nLoc(x2))
  expect_equal(nInd(x2), nInd(x))
  expect_true(all(ploidy(x2) == 2))
  expect_length(x2@other$history, length(x@other$history) + 1)
})

test_that("gl.filter.callrate method='ind' removes individuals below threshold (SNP)", {
  x <- testset.gl
  icr <- 1 - rowSums(is.na(as.matrix(x))) / nLoc(x)
  out <- capture.output(x2 <- gl.filter.callrate(x, method = "ind", threshold = 0.9,
                                                 plot.display = FALSE, verbose = 0))
  expect_equal(nInd(x2), sum(icr >= 0.9))
  expect_equal(nLoc(x2), nLoc(x))
  expect_true(all(indNames(x2) %in% indNames(x)[icr >= 0.9]))
})

test_that("gl.filter.callrate method='pop' retains only loci passing in every population", {
  x <- testset.gl
  out <- capture.output(x2 <- gl.filter.callrate(x, method = "pop", threshold = 0.8,
                                                 plot.display = FALSE, verbose = 0))
  expect_true(nLoc(x2) <= nLoc(x))
  expect_equal(nInd(x2), nInd(x))
  expect_equal(nrow(x2@other$loc.metrics), nLoc(x2))
  # independent check: every retained locus passes the threshold in every population
  x_rc <- x2
  for (p in levels(pop(x_rc))) {
    sub <- x_rc[pop(x_rc) == p, ]
    cr <- 1 - colSums(is.na(as.matrix(sub))) / nInd(sub)
    expect_true(all(cr >= 0.8), label = paste("population", p))
  }
})

test_that("gl.filter.callrate works on SilicoDArT (loc and ind)", {
  x <- testset.gs
  out <- capture.output(x2 <- gl.filter.callrate(x, method = "loc", threshold = 0.9,
                                                 plot.display = FALSE, verbose = 0))
  expect_true(all(ploidy(x2) == 1))
  out2 <- capture.output(x3 <- gl.filter.callrate(x, method = "ind", threshold = 0.8,
                                                  plot.display = FALSE, verbose = 0))
  expect_true(all(ploidy(x3) == 1))
})

test_that("monomorphs flag NULL no longer crashes (F2 fix)", {
  x <- testset.gl
  x_null <- x
  x_null@other$loc.metrics.flags$monomorphs <- NULL
  expect_error(
    capture.output(gl.filter.callrate(x_null, method = "loc", threshold = 0.9,
                                      plot.display = FALSE, verbose = 0)),
    NA
  )
})

test_that("pop branch preserves locus metrics intact (F1 fix)", {
  x <- testset.gl
  out <- capture.output(xp <- gl.filter.callrate(x, method = "pop", threshold = 0.8,
                                                 plot.display = FALSE, verbose = 0))
  expect_equal(sum(!complete.cases(xp@other$loc.metrics)), 0)
  # metrics rows must be the source rows for exactly the retained loci
  keep <- locNames(x) %in% locNames(xp)
  expect_equal(as.character(xp@other$loc.metrics$AlleleID),
               as.character(x@other$loc.metrics$AlleleID[keep]))
})

test_that("at-threshold individuals retained and NOT listed as deleted (F3 fix)", {
  x <- testset.gl
  icr <- 1 - rowSums(is.na(as.matrix(x))) / nLoc(x)
  thr <- sort(unique(icr))[3]
  out2 <- capture.output(x2 <- gl.filter.callrate(x, method = "ind", threshold = thr,
                                                  plot.display = FALSE, verbose = 2))
  at_thr <- indNames(x)[icr == thr]
  expect_gt(length(at_thr), 0)
  del_block <- paste(out2, collapse = " ")
  expect_true(all(at_thr %in% indNames(x2)))
  expect_false(any(sapply(at_thr, function(n) grepl(n, del_block, fixed = TRUE))))
})

test_that("verbose = 0 is fully silent even when individuals are deleted (F5 fix)", {
  x <- testset.gl
  out3 <- capture.output(gl.filter.callrate(x, method = "ind", threshold = 0.9,
                                            plot.display = FALSE, verbose = 0))
  expect_length(out3, 0)
})

test_that("monomorphs flag invalidated after ind-filtering with mono.rm=FALSE (F4 fix)", {
  x <- testset.gl
  x_mono <- gl.filter.monomorphs(x, verbose = 0)
  expect_true(x_mono@other$loc.metrics.flags$monomorphs)
  out4 <- capture.output(x4 <- gl.filter.callrate(x_mono, method = "ind", threshold = 0.9,
                                                  plot.display = FALSE, verbose = 0))
  expect_false(x4@other$loc.metrics.flags$monomorphs)
})
