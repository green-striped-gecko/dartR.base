# Characterization tests for gl.filter.ld
# Baseline snapshotted before review (review-gl.filter.ld).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

prep_platypus_mapped <- function() {
  x <- platypus.gl
  x <- gl.filter.callrate(x, threshold = 1, verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  x$position <- x$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
  x$chromosome <- as.factor(x$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1)
  x
}

ld_report_platypus <- function(x) {
  gl.report.ld.map(x, ld.max.pairwise = 10000000,
                   plot.display = FALSE, verbose = 0)
}

test_that("gl.filter.ld baseline on the platypus LD report", {
  x <- prep_platypus_mapped()
  res <- ld_report_platypus(x)
  out <- capture.output(f1 <- gl.filter.ld(x, ld.report = res, verbose = 0))
  expect_length(out, 0)
  expect_equal(nLoc(f1), 380)
  expect_equal(nrow(f1@other$loc.metrics), nLoc(f1))
  f2 <- gl.filter.ld(x, ld.report = res, pop.limit = 1, verbose = 0)
  expect_equal(nLoc(f2), 349)
  expect_equal(nrow(f2@other$loc.metrics), nLoc(f2))
})

test_that("gl.filter.ld default pop.limit is half the populations in the report", {
  # [approved diff F6] the default is now computed explicitly from
  # ld.report$pop (previously the same value arose from lazy evaluation of
  # ceiling(nPop(x)/2) after x had been subset to the report populations).
  # testset.gl has 30 pops; only 3 clear ind.limit = 10 -> default is 2.
  tt <- gl.filter.allna(testset.gl, verbose = 0)
  tt <- gl.filter.monomorphs(tt, verbose = 0)
  tt <- tt[, 1:80]
  tt@other$loc.metrics <- tt@other$loc.metrics[1:80, ]
  out <- capture.output(
    rt <- gl.report.ld.map(tt, plot.display = FALSE, verbose = 0)
  )
  expect_equal(length(unique(rt$pop)), 3)
  fdef <- gl.filter.ld(tt, ld.report = rt, threshold = 0.2, verbose = 0)
  f2 <- gl.filter.ld(tt, ld.report = rt, threshold = 0.2, pop.limit = 2,
                     verbose = 0)
  expect_equal(nLoc(fdef), nLoc(f2))
  expect_equal(nLoc(fdef), 79)
})

test_that("gl.filter.ld drops the partner of an already-dropped locus", {
  # [pinned behaviour — change 1 REJECTED 2026-09-10] comparisons are
  # strictly pairwise and sequential: L3 is removed although its only LD
  # partner L2 was already removed. This is now documented in @details.
  x <- prep_platypus_mapped()
  L <- locNames(x)[1:3]
  rep_chain <- data.frame(
    pop = "TENTERFIELD",
    chr = "c1", pos_loc_a = c(1, 2), pos_loc_b = c(2, 3),
    ld.stat = c(0.9, 0.9), distance = c(1, 1),
    locus_a.snp.name = c(L[1], L[2]),
    locus_a.stat.keep = c(3, 2),
    locus_b.snp.name = c(L[2], L[3]),
    locus_b.stat.keep = c(2, 1),
    locus_a_b = c(paste0(L[1], "_", L[2]), paste0(L[2], "_", L[3]))
  )
  fc <- gl.filter.ld(x, ld.report = rep_chain, pop.limit = 1, verbose = 0)
  dropped <- setdiff(locNames(x), locNames(fc))
  expect_true(L[2] %in% dropped)
  expect_true(L[3] %in% dropped)
  expect_false(L[1] %in% dropped)
})

test_that("gl.filter.ld no-pairs branch is silent at verbose = 0", {
  # [approved diff F2] the "No pair of loci ..." message gates at verbose >= 1
  x <- prep_platypus_mapped()
  res <- ld_report_platypus(x)
  out <- capture.output(
    f0 <- gl.filter.ld(x, ld.report = res, threshold = 1.5, verbose = 0)
  )
  expect_length(out, 0)
  expect_equal(nLoc(f0), nLoc(x))
  out1 <- capture.output(
    f0 <- gl.filter.ld(x, ld.report = res, threshold = 1.5, verbose = 1)
  )
  expect_true(any(grepl("No pair of loci", out1)))
})

test_that("gl.filter.ld appends a single history entry per call", {
  # [approved diff F3] previously two entries: one from the internal
  # gl.drop.loc call (exposing internal variable names) and one of its own
  x <- prep_platypus_mapped()
  res <- ld_report_platypus(x)
  h <- length(x@other$history)
  f1 <- gl.filter.ld(x, ld.report = res, verbose = 0)
  expect_equal(length(f1@other$history), h + 1)
  expect_true(grepl("gl.filter.ld",
                    deparse(f1@other$history[[h + 1]])[1]))
})

test_that("gl.filter.ld tolerates loc.metrics.flags$monomorphs absent", {
  # [approved diff F4] objects not built by dartR (flag absent) previously
  # crashed with "argument is of length zero"; now they take the warning
  # path and filter normally
  x <- prep_platypus_mapped()
  res <- ld_report_platypus(x)
  z <- x
  z@other$loc.metrics.flags$monomorphs <- NULL
  expect_no_error(f <- gl.filter.ld(z, ld.report = res, verbose = 0))
  expect_equal(nLoc(f), 380)
})

test_that("gl.filter.ld rejects an ld.report without the required columns", {
  # [approved diff F5] fail-fast validation instead of obscure downstream
  # errors
  x <- prep_platypus_mapped()
  expect_error(
    gl.filter.ld(x, ld.report = data.frame(a = 1), verbose = 0),
    "gl.report.ld.map"
  )
})

test_that("gl.filter.ld returns invisibly", {
  x <- prep_platypus_mapped()
  res <- ld_report_platypus(x)
  v <- withVisible(gl.filter.ld(x, ld.report = res, verbose = 0))
  expect_false(v$visible)
})
