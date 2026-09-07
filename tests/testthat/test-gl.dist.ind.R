# Characterization tests for gl.dist.ind
# Baseline snapshotted before review (population-distance chain review,
# dev at ddaed27). Pins current behaviour, defects included; assertions
# tagged [approved Fn] were flipped when the matching approved finding
# was applied on branch review-gl.dist.ind.

# utils.dist.ind.snp (fixed on the open PR #315, not yet on dev) leaks
# one progress line at verbose = 0 on the SNP paths. Filter it so the
# silence assertions test gl.dist.ind's own gating; the filter becomes a
# no-op once #315 merges.
strip_engine_leak <- function(o) {
  o[!grepl("Calculating the (un)?scaled distance matrix|Calculating the distance matrix", o)]
}

test_that("SNP wrapper output equals the utils.dist.ind.snp engine", {
  gli <- testset.gl[1:8, 1:60]
  for (m in c("euclidean", "simple", "manhattan", "czekanowski",
              "absolute")) {
    Dw <- gl.dist.ind(gli, method = m, plot.display = FALSE, verbose = 0)
    De <- as.dist(utils.dist.ind.snp(gli, method = m, verbose = 0))
    expect_equal(as.vector(Dw), as.vector(De), info = m)
  }
  Ds <- gl.dist.ind(gli, method = "euclidean", scale = TRUE,
                    plot.display = FALSE, verbose = 0)
  De <- as.dist(utils.dist.ind.snp(gli, method = "euclidean",
                                   scale = TRUE, verbose = 0))
  expect_equal(as.vector(Ds), as.vector(De))
})

test_that("SNP anchors on testset.gl[1:8, 1:60]", {
  gli <- testset.gl[1:8, 1:60]
  De <- gl.dist.ind(gli, method = "euclidean", plot.display = FALSE,
                    verbose = 0)
  expect_equal(round(as.vector(De)[1:5], 6),
               c(1.414214, 0, 0, 0, 1))
  Dm <- gl.dist.ind(gli, method = "manhattan", plot.display = FALSE,
                    verbose = 0)
  expect_equal(round(as.vector(Dm)[1:5], 6),
               c(0.018519, 0, 0, 0, 0.009434))
})

test_that("scale=TRUE is ignored for non-euclidean methods", {
  # [approved F3, docs-only] behaviour unchanged; the docs now state
  # that scale applies to euclidean only and default FALSE
  gli <- testset.gl[1:8, 1:60]
  D1 <- gl.dist.ind(gli, method = "simple", scale = TRUE,
                    plot.display = FALSE, verbose = 0)
  D0 <- gl.dist.ind(gli, method = "simple", scale = FALSE,
                    plot.display = FALSE, verbose = 0)
  expect_equal(as.vector(D1), as.vector(D0))
})

test_that("silico wrapper output equals the utils.dist.binary engine", {
  gsi <- testset.gs[1:8, 1:60]
  for (m in c("euclidean", "simple", "jaccard")) {
    Dw <- gl.dist.ind(gsi, method = m, plot.display = FALSE, verbose = 0)
    De <- as.dist(utils.dist.binary(gsi, method = m, verbose = 0))
    expect_equal(as.vector(Dw), as.vector(De), info = m)
  }
  Ds <- gl.dist.ind(gsi, method = "simple", plot.display = FALSE,
                    verbose = 0)
  expect_equal(round(as.vector(Ds)[1:5], 6),
               c(0.089286, 0.104167, 0, 0.148148, 0.218182))
  Dj <- gl.dist.ind(gsi, method = "jaccard", plot.display = FALSE,
                    verbose = 0)
  expect_equal(round(as.vector(Dj)[1:5], 6),
               c(0.227273, 0.333333, 0, 0.347826, 0.5))
})

test_that("'sorensen' is accepted and returns the engine's Sorensen", {
  # [approved F1] sorensen added to the validation list; previously it
  # was coerced to simple matching with an ungated warning
  gsi <- testset.gs[1:8, 1:60]
  o <- capture.output(Dso <- gl.dist.ind(gsi, method = "sorensen",
                                         plot.display = FALSE,
                                         verbose = 0))
  expect_equal(length(o), 0)
  Deng <- as.dist(utils.dist.binary(gsi, method = "sorensen", verbose = 0))
  expect_equal(as.vector(Dso), as.vector(Deng))
  Dsi <- gl.dist.ind(gsi, method = "simple", plot.display = FALSE,
                     verbose = 0)
  expect_false(isTRUE(all.equal(as.vector(Dso), as.vector(Dsi))))
})

test_that("an unknown SNP method falls back quietly at verbose 0", {
  # [approved F2] the fallback warnings are now gated at verbose >= 1;
  # the fallback to euclidean itself is unchanged
  gli <- testset.gl[1:8, 1:60]
  o <- capture.output(Dx <- gl.dist.ind(gli, method = "foo",
                                        plot.display = FALSE,
                                        verbose = 0))
  expect_equal(length(strip_engine_leak(o)), 0)
  De <- gl.dist.ind(gli, method = "euclidean", plot.display = FALSE,
                    verbose = 0)
  expect_equal(as.vector(Dx), as.vector(De))
  # at verbose = 1 the warning prints
  o1 <- capture.output(invisible(gl.dist.ind(gli, method = "foo",
                                             plot.display = FALSE,
                                             verbose = 1)))
  expect_true(any(grepl("not in the list of options", o1)))
})

test_that("type='matrix' returns a full symmetric named matrix", {
  gli <- testset.gl[1:8, 1:60]
  Dm <- gl.dist.ind(gli, method = "euclidean", type = "matrix",
                    plot.display = FALSE, verbose = 0)
  expect_true(is.matrix(Dm))
  expect_true(isSymmetric(unname(Dm)))
  expect_true(all(diag(Dm) == 0))
  expect_identical(rownames(Dm), indNames(gli))
  # [approved F5] type is normalised with tolower() and validated:
  # 'Matrix' now returns a matrix; an unrecognised type stops
  Dt <- gl.dist.ind(gli, method = "euclidean", type = "Matrix",
                    plot.display = FALSE, verbose = 0)
  expect_true(is.matrix(Dt))
  expect_error(gl.dist.ind(gli, method = "euclidean", type = "foo",
                           plot.display = FALSE, verbose = 0),
               "dist")
})

test_that("verbose=0 is silent and plot.file works with display off", {
  gli <- testset.gl[1:8, 1:60]
  o <- capture.output(invisible(gl.dist.ind(gli, method = "euclidean",
                                            verbose = 0)))
  expect_equal(length(strip_engine_leak(o)), 0)
  expect_no_error(gl.dist.ind(gli, method = "euclidean",
                              plot.display = FALSE, plot.file = "tstI",
                              plot.dir = tempdir(), verbose = 0))
})

test_that("an all-missing individual yields NA distances with a warning", {
  gna <- testset.gl[1:6, 1:60]
  mna <- as.matrix(gna)
  mna[3, ] <- NA
  g2 <- new("genlight", gen = mna, ploidy = 2)
  indNames(g2) <- indNames(gna)
  locNames(g2) <- locNames(gna)
  pop(g2) <- pop(gna)
  g2 <- gl.compliance.check(g2, verbose = 0)
  D <- gl.dist.ind(g2, method = "euclidean", plot.display = FALSE,
                   verbose = 0)
  M <- as.matrix(D)
  expect_true(all(is.na(M[3, -3])))
  # [approved F6] NA distances are counted and warned at verbose >= 1
  # (still silent at verbose = 0)
  o1 <- capture.output(invisible(
    gl.dist.ind(g2, method = "euclidean", plot.display = FALSE,
                verbose = 1)))
  expect_true(any(grepl("distances are NA", o1)))
})
