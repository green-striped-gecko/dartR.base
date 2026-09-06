# Characterization tests for gl.tree.nj
# Baseline snapshotted before review (review-gl.tree.nj) against
# upstream/dev ddaed27, then updated in Phase C: assertions marked
# [approved Fn] pin the repaired behaviour approved 2026-09-06 (see
# function-review/reports/dartR.base/gl.tree.nj.md).

# Route all graphics to a null device per call so tests stay headless.
quiet_tree <- function(...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  gl.tree.nj(...)
}

# Independent reference: replicate the documented computation (Euclidean
# distance on per-population allele frequencies) directly from the
# genotype matrix. Frequencies exclude missing genotypes and are scaled
# by ploidy [approved F1, F7]; 4 dp rounding as in the function.
ref_dist <- function(gl, na.rm = TRUE, pl = 2) {
  m <- as.matrix(gl)
  f <- apply(m, 2, tapply, pop(gl), function(e) mean(e, na.rm = na.rm) / pl)
  as.dist(round(as.matrix(dist(f)), 4))
}

test_that("returns phylo over populations; matches independent recomputation", {
  xcopy <- testset.gl
  o <- capture.output(v <- withVisible(quiet_tree(testset.gl, verbose = 0)))
  expect_length(o, 0)                      # silent at verbose = 0 (text side)
  expect_true(v$visible)                   # current behaviour: visible return
  tr <- v$value
  expect_s3_class(tr, "phylo")
  expect_identical(xcopy, testset.gl)      # input object untouched
  expect_setequal(tr$tip.label, popNames(testset.gl))
  expect_length(tr$tip.label, 30L)
  # Independent verification: same tree as ape::nj on the recomputed
  # Euclidean frequency distance (RF distance 0, identical edge sums)
  tref <- ape::nj(ref_dist(testset.gl))
  expect_equal(ape::dist.topo(ape::unroot(tr), ape::unroot(tref))[1], 0)
  expect_equal(sum(tr$edge.length), sum(tref$edge.length), tolerance = 1e-8)
  expect_equal(sum(tr$edge.length), 18.23083, tolerance = 1e-4) # dartR.data 1.2.5
  # [approved F1] frequencies now exclude missing genotypes (na.rm = TRUE,
  # matching gl.dist.pop); the old no-na.rm matrix gives a materially
  # different topology
  told <- ape::nj(ref_dist(testset.gl, na.rm = FALSE))
  expect_equal(ape::dist.topo(ape::unroot(tr), ape::unroot(told))[1], 32)
})

test_that("user-supplied dist.matrix is used as-is", {
  d <- ref_dist(testset.gl)
  tr <- quiet_tree(testset.gl, dist.matrix = d, verbose = 0)
  expect_equal(ape::dist.topo(ape::unroot(tr), ape::unroot(ape::nj(d)))[1], 0)
  # [approved F8] a plain matrix is coerced with as.dist(), so the upgma
  # path no longer fails opaquely inside hclust
  tm <- quiet_tree(testset.gl, dist.matrix = as.matrix(d), method = "upgma",
                   verbose = 0)
  href <- ape::as.phylo(stats::hclust(d, method = "average"))
  expect_equal(ape::dist.topo(ape::unroot(tm), ape::unroot(href))[1], 0)
})

test_that("method: UPGMA accepted; ugpma retained as synonym [approved F3]", {
  trnj <- quiet_tree(testset.gl, verbose = 0)
  tru <- quiet_tree(testset.gl, method = "UPGMA", verbose = 0)
  expect_s3_class(tru, "phylo")
  href <- ape::as.phylo(stats::hclust(ref_dist(testset.gl), method = "average"))
  expect_equal(ape::dist.topo(ape::unroot(tru), ape::unroot(href))[1], 0)
  expect_gt(ape::dist.topo(ape::unroot(tru), ape::unroot(trnj))[1], 0)
  # legacy misspelling produces the identical UPGMA tree
  trl <- quiet_tree(testset.gl, method = "ugpma", verbose = 0)
  expect_equal(ape::dist.topo(ape::unroot(trl), ape::unroot(tru))[1], 0)
  # unknown method still falls back to nj, but the warning is gated:
  # silent at verbose = 0, printed at verbose >= 1
  o0 <- capture.output(trw <- quiet_tree(testset.gl, method = "banana",
                                         verbose = 0))
  expect_length(o0, 0)
  expect_equal(ape::dist.topo(ape::unroot(trw), ape::unroot(trnj))[1], 0)
  o1 <- capture.output(invisible(quiet_tree(testset.gl, method = "banana",
                                            verbose = 1)))
  expect_match(paste(o1, collapse = " "), "method must be one of 'nj' or 'upgma'")
})

test_that("by.pop = FALSE returns an individual-labelled tree [approved F2]", {
  tr <- quiet_tree(testset.gl, by.pop = FALSE, verbose = 0)
  expect_s3_class(tr, "phylo")
  expect_setequal(tr$tip.label, indNames(testset.gl))
  expect_length(tr$tip.label, nInd(testset.gl))
})

test_that("as.pop relabels tips from an ind.metrics column", {
  tr <- quiet_tree(testset.gl, as.pop = "sex", verbose = 0)
  expect_setequal(tr$tip.label, c("Female", "Male", "Unknown"))
  # [approved F6] error message now cites ind.metrics for the
  # ind.metrics lookup
  expect_error(quiet_tree(testset.gl, as.pop = "nosuchcol", verbose = 0),
               "ind.metrics")
})

test_that("outgroup re-roots; treefile written; bad outgroup errors", {
  tf <- tempfile(fileext = ".tre")
  tr <- quiet_tree(testset.gl, outgroup = "EmmacBrisWive", treefile = tf,
                   verbose = 0)
  expect_true(file.exists(tf))
  expect_s3_class(ape::read.tree(tf), "phylo")
  # [approved F5] rooting now resolves the root, so the returned tree is
  # rooted in ape's sense
  expect_true(ape::is.rooted(tr))
  expect_error(quiet_tree(testset.gl, outgroup = "NOSUCHPOP", verbose = 0),
               "outgroup")
  unlink(tf)
})

test_that("SilicoDArT frequencies scaled by ploidy 1 [approved F7]", {
  tr <- quiet_tree(testset.gs, verbose = 0)
  expect_s3_class(tr, "phylo")
  expect_setequal(tr$tip.label, popNames(testset.gs))
  tref <- ape::nj(ref_dist(testset.gs, pl = 1))
  expect_equal(ape::dist.topo(ape::unroot(tr), ape::unroot(tref))[1], 0)
  expect_equal(sum(tr$edge.length), sum(tref$edge.length), tolerance = 1e-8)
  expect_equal(sum(tr$edge.length), 47.33798, tolerance = 1e-4) # dartR.data 1.2.5
})

test_that("population-count edge cases", {
  p <- popNames(testset.gl)
  g2 <- gl.keep.pop(testset.gl, pop.list = p[1:2], verbose = 0)
  expect_error(quiet_tree(g2, verbose = 0),
               "less than 3 observations")     # raw ape error, no dartR guard
  g3 <- gl.keep.pop(testset.gl, pop.list = p[1:3], verbose = 0)
  t3 <- quiet_tree(g3, verbose = 0)
  expect_length(t3$tip.label, 3L)
})

test_that("invalid type fails fast; plot failure keeps the tree [approved F4]", {
  # type is validated during error checking, before any computation
  expect_error(quiet_tree(testset.gl, type = "banana", verbose = 0),
               "phylogram")
  # a failure inside the plot call no longer destroys the computed tree
  local_mocked_bindings(plot.phylo = function(...) stop("no device"),
                        .package = "ape")
  o <- capture.output(tr <- quiet_tree(testset.gl, verbose = 1))
  expect_s3_class(tr, "phylo")
  expect_match(paste(o, collapse = " "), "could not be plotted")
  # at verbose = 0 the plot is gated off entirely (plot.display forced
  # FALSE), so the mocked failure is never reached
  o0 <- capture.output(tr0 <- quiet_tree(testset.gl, verbose = 0))
  expect_length(o0, 0)
  expect_s3_class(tr0, "phylo")
})
