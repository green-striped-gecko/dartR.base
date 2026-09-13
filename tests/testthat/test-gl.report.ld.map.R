# Characterization tests for gl.report.ld.map
# Baseline snapshotted before review (review-gl.report.ld.map).
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

prep_testset_unmapped <- function() {
  tt <- gl.filter.allna(testset.gl, verbose = 0)
  tt <- gl.filter.monomorphs(tt, verbose = 0)
  tt <- tt[, 1:80]
  tt@other$loc.metrics <- tt@other$loc.metrics[1:80, ]
  tt
}

test_that("gl.report.ld.map mapped baseline on platypus.gl", {
  x <- prep_platypus_mapped()
  expect_equal(nLoc(x), 383)
  res <- gl.report.ld.map(x, ld.max.pairwise = 10000000,
                          plot.display = FALSE, verbose = 0)
  expect_equal(nrow(res), 496)
  expect_equal(colnames(res), c(
    "pop", "chr", "pos_loc_a", "pos_loc_b", "ld.stat", "distance",
    "locus_a.snp.name", "locus_a.stat.keep", "locus_b.snp.name",
    "locus_b.stat.keep", "locus_a_b"
  ))
  expect_equal(as.vector(table(res$pop)), c(141, 173, 182))
  expect_equal(max(res$ld.stat), 1, tolerance = 1e-8)
  expect_equal(mean(res$ld.stat), 0.0775235, tolerance = 1e-5)
  # input untouched (report family)
  x2 <- prep_platypus_mapped()
  expect_identical(x, x2)
})

test_that("gl.report.ld.map reported pairs carry correct labels and positions", {
  # independent recomputation: snpStats R.squared on the largest population,
  # keyed by the reported locus names, reproduces ld.stat exactly
  x <- prep_platypus_mapped()
  res <- gl.report.ld.map(x, ld.max.pairwise = 10000000,
                          plot.display = FALSE, verbose = 0)
  pops <- seppop(x)
  big <- pops[[which.max(sapply(pops, nInd))]]
  big <- gl.filter.maf(big, threshold = 0.05, verbose = 0)
  m <- as.matrix(big)
  r <- m + 1; r[is.na(r)] <- 0
  sm <- methods::new("SnpMatrix",
                     matrix(as.raw(r), nrow = nrow(m),
                            dimnames = list(rownames(m), colnames(m))))
  suppressWarnings(r2m <- as.matrix(snpStats::ld(sm, sm, stats = "R.squared")))
  sub <- res[res$pop == popNames(big), ]
  set.seed(42)
  take <- sub[sample(nrow(sub), 8), ]
  for (i in seq_len(nrow(take))) {
    a <- as.character(take$locus_a.snp.name[i])
    b <- as.character(take$locus_b.snp.name[i])
    ind <- max(r2m[a, b], r2m[b, a], na.rm = TRUE)
    expect_equal(take$ld.stat[i], ind, tolerance = 1e-9)
  }
  posmap <- stats::setNames(x$position, locNames(x))
  expect_true(all(posmap[as.character(take$locus_a.snp.name)] == take$pos_loc_a))
  expect_true(all(posmap[as.character(take$locus_b.snp.name)] == take$pos_loc_b))
})

test_that("gl.report.ld.map unmapped auto-detect on testset.gl subset", {
  tt <- prep_testset_unmapped()
  out <- capture.output(
    res <- gl.report.ld.map(tt, plot.display = FALSE, verbose = 0)
  )
  expect_equal(nrow(res), 17)
  expect_setequal(unique(as.character(res$pop)),
                  c("EmmacBurnBara", "EmmacMaclGeor", "EmmacMDBForb"))
  # explicit ld.max.pairwise = NULL takes the same path
  out2 <- capture.output(
    res2 <- gl.report.ld.map(tt, ld.max.pairwise = NULL,
                             plot.display = FALSE, verbose = 0)
  )
  expect_equal(res$ld.stat, res2$ld.stat)
  # [approved diff F3] fully silent at verbose = 0
  expect_length(out, 0)
  # [approved diff F3] skipped-population warnings at verbose >= 1, the
  # no-chromosome note at verbose >= 2
  out1 <- capture.output(
    res3 <- gl.report.ld.map(tt, plot.display = FALSE, verbose = 1)
  )
  expect_true(any(grepl("Skipping population", out1)))
  expect_false(any(grepl("chromosome/position", out1)))
  out3 <- capture.output(
    res4 <- gl.report.ld.map(tt, plot.display = FALSE, verbose = 2)
  )
  expect_true(any(grepl("chromosome/position", out3)))
})

test_that("gl.report.ld.map with a signed statistic returns no negative values", {
  # [pinned truncation — change 1 REJECTED 2026-09-10] pairs with a zero or
  # negative statistic are not reported; this is now documented behaviour
  x <- prep_platypus_mapped()
  res_r <- gl.report.ld.map(x, ld.max.pairwise = 10000000, ld.stat = "R",
                            plot.display = FALSE, verbose = 0)
  expect_equal(nrow(res_r), 234)
  expect_true(min(res_r$ld.stat) > 0)
})

test_that("gl.report.ld.map saves the plot with plot.display = FALSE", {
  # [approved diff F2] plots are built whenever displayed OR saved; this
  # call used to fail with "object 'p4' not found"
  x <- prep_platypus_mapped()
  pd <- file.path(tempdir(), "ldmap-plots")
  dir.create(pd, showWarnings = FALSE)
  expect_no_error(
    res <- gl.report.ld.map(x, ld.max.pairwise = 10000000,
                            plot.display = FALSE, plot.file = "ldtest",
                            plot.dir = pd, verbose = 0)
  )
  expect_true(length(list.files(pd, pattern = "ldtest")) > 0)
})

test_that("gl.report.ld.map rejects SilicoDArT", {
  # [approved diff F6] rejected at the entry datatype check (accept = "SNP")
  gs <- testset.gs[, 1:10]
  gs@other$loc.metrics <- testset.gs@other$loc.metrics[1:10, ]
  expect_error(
    gl.report.ld.map(gs, plot.display = FALSE, verbose = 0),
    "SilicoDArT"
  )
})
