# Characterization tests for utils.recalc.avgpic
# Baseline snapshotted before review (recalc-battery review, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.

test_that("AvgPIC/PIC components match an independent recomputation", {
  sub <- gl.drop.pop(testset.gl, pop.list = popNames(testset.gl)[1:10],
                     verbose = 0)
  gm <- as.matrix(sub)
  o <- capture.output(r <- utils.recalc.avgpic(sub, verbose = 0))
  expect_equal(length(o), 0)
  c0 <- colSums(gm == 0, na.rm = TRUE); c1 <- colSums(gm == 1, na.rm = TRUE)
  c2 <- colSums(gm == 2, na.rm = TRUE); ct <- c0 + c1 + c2
  orr <- (c0 + c1) / ct; ors <- (c1 + c2) / ct
  picr <- 1 - (orr^2 + (1 - orr)^2); pics <- 1 - (ors^2 + (1 - ors)^2)
  expect_equal(unname(r@other$loc.metrics$AvgPIC), unname((picr + pics) / 2))
  expect_equal(unname(r@other$loc.metrics$OneRatioRef), unname(orr))
})

test_that("SilicoDArT branch computes OneRatio and PIC", {
  o <- capture.output(r <- utils.recalc.avgpic(testset.gs, verbose = 0))
  onr <- colMeans(as.matrix(testset.gs) == 1, na.rm = TRUE)
  expect_equal(unname(r@other$loc.metrics$OneRatio), unname(onr))
  expect_equal(unname(r@other$loc.metrics$PIC),
               unname(1 - (onr^2 + (1 - onr)^2)))
})

test_that("an object without a flags slot is handled", {
  sub <- gl.drop.pop(testset.gl, pop.list = popNames(testset.gl)[1:10],
                     verbose = 0)
  sub@other$loc.metrics.flags <- NULL
  # [approved diff B2] baseline: the monomorphs check used `flag == FALSE`,
  # which is `if (logical(0))` when the slot is absent.
  expect_error(capture.output(utils.recalc.avgpic(sub, verbose = 0)),
    NA)  # [approved diff B2]
})

