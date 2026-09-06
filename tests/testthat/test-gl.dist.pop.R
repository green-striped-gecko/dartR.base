# Characterization tests for gl.dist.pop
# Baseline snapshotted before review (population-distance chain review,
# dev at ddaed27). Pins current behaviour, defects included; assertions
# tagged [approved Fn] were flipped when the matching approved finding
# was applied on branch review-gl.dist.pop. Anchor values were computed
# on this machine at the reviewed state and are 6 dp rounded.

glp4 <- function() gl.keep.pop(testset.gl, popNames(testset.gl)[1:4],
                               verbose = 0)
gsp4 <- function() gl.keep.pop(testset.gs, popNames(testset.gs)[1:4],
                               verbose = 0)

test_that("SNP method anchors on 4 pops of testset.gl", {
  glp <- glp4()
  anchors <- list(
    euclidean  = c(1.586012, 1.745165, 1.286407, 0.945351, 1.166444,
                   1.504021),
    nei        = c(0.010975, 0.012982, 0.007362, 0.003887, 0.006184,
                   0.010153),
    reynolds   = c(1.217084, 1.080242, 1.019249, 0.635415, 0.941806,
                   1.026265),
    chord      = c(0.096749, 0.096955, 0.081307, 0.068248, 0.082704,
                   0.095207),
    `fixed-diff` = c(0.04, 0, 0, 0, 0.04, 0)
  )
  for (m in names(anchors)) {
    D <- suppressWarnings(gl.dist.pop(glp, method = m,
                                      plot.display = FALSE, verbose = 0))
    expect_s3_class(D, "dist")
    expect_equal(round(as.vector(D), 6), anchors[[m]], info = m)
    expect_equal(attr(D, "Size"), 4, info = m)
  }
  Ds <- gl.dist.pop(glp, method = "euclidean", scale = TRUE,
                    plot.display = FALSE, verbose = 0)
  expect_equal(round(as.vector(Ds), 6),
               c(0.051952, 0.056325, 0.042597, 0.030966, 0.039055,
                 0.049913))
})

test_that("frequency methods recompute exactly from gl.allele.freq cells", {
  glp <- glp4()
  m <- as.matrix(glp)
  P <- t(sapply(popNames(glp), function(p) {
    round(colMeans(m[as.character(pop(glp)) == p, , drop = FALSE],
                   na.rm = TRUE) / 2 * 100, 2) / 100
  }))
  D <- as.matrix(gl.dist.pop(glp, method = "euclidean",
                             plot.display = FALSE, verbose = 0))
  for (i in 1:3) for (j in (i + 1):4) {
    sq <- (P[i, ] - P[j, ])^2
    sq <- sq[!is.na(sq)]
    expect_equal(D[j, i], sqrt(sum(sq)))
  }
})

test_that("fixed-diff equals gl.fixed.diff $pcfd / 100", {
  glp <- glp4()
  D <- suppressWarnings(gl.dist.pop(glp, method = "fixed-diff",
                                    plot.display = FALSE, verbose = 0))
  fd <- suppressWarnings(gl.fixed.diff(glp, verbose = 0))[[3]]
  expect_equal(as.vector(D), as.vector(fd) / 100)
})

test_that("distances follow their labels for any pop level order", {
  # [approved F1] the frequency matrix rows are re-anchored to
  # popNames(x) by name after the dcast, so reversed or shuffled
  # population factor levels return the same distances under the same
  # labels as the alphabetical-level object.
  glp <- glp4()
  D0 <- as.matrix(gl.dist.pop(glp, method = "euclidean",
                              plot.display = FALSE, verbose = 0))
  # reversed levels
  glr <- glp
  pop(glr) <- factor(as.character(pop(glp)),
                     levels = rev(sort(unique(as.character(pop(glp))))))
  Dr <- as.matrix(gl.dist.pop(glr, method = "euclidean",
                              plot.display = FALSE, verbose = 0))
  expect_true(isTRUE(all.equal(Dr[rownames(D0), colnames(D0)], D0,
                               tolerance = 1e-9)))
  # shuffled levels
  gls <- glp
  set.seed(11)
  pop(gls) <- factor(as.character(pop(glp)),
                     levels = sample(unique(as.character(pop(glp)))))
  Ds <- as.matrix(gl.dist.pop(gls, method = "nei",
                              plot.display = FALSE, verbose = 0))
  Dn <- as.matrix(gl.dist.pop(glp, method = "nei",
                              plot.display = FALSE, verbose = 0))
  expect_true(isTRUE(all.equal(Ds[rownames(Dn), colnames(Dn)], Dn,
                               tolerance = 1e-9)))
})

test_that("an unknown method stops with a fatal error", {
  # [approved F4] previously an ungated non-fatal "Fatal Error" printed
  # at verbose 0 and the function silently ran euclidean
  glp <- glp4()
  expect_error(gl.dist.pop(glp, method = "neii", plot.display = FALSE,
                           verbose = 0),
               "not among those")
})

test_that("type='matrix' returns a full symmetric zero-diagonal matrix", {
  # [approved F5] previously the SNP frequency methods returned the
  # lower triangle only, with the upper triangle and diagonal NA
  glp <- glp4()
  Dm <- gl.dist.pop(glp, method = "euclidean", type = "matrix",
                    plot.display = FALSE, verbose = 0)
  expect_true(isSymmetric(unname(Dm)))
  expect_true(all(diag(Dm) == 0))
  expect_false(anyNA(Dm))
  Dd <- as.matrix(gl.dist.pop(glp, method = "euclidean", type = "dist",
                              plot.display = FALSE, verbose = 0))
  expect_equal(unname(Dm), unname(Dd))
})

test_that("plot.file with plot.display=FALSE saves and returns the dist", {
  # [approved F3] the plot objects are now built whenever displayed OR
  # saved; previously this call died with "object 'p3' not found" after
  # the distances were computed
  glp <- glp4()
  D <- gl.dist.pop(glp, method = "euclidean", plot.display = FALSE,
                   plot.file = "tstX", plot.dir = tempdir(), verbose = 0)
  expect_s3_class(D, "dist")
  expect_true(file.exists(file.path(tempdir(), "tstX.RDS")))
})

test_that("verbose=0 is silent for valid methods", {
  glp <- glp4()
  o <- capture.output(invisible(gl.dist.pop(glp, method = "nei",
                                            verbose = 0)))
  expect_equal(length(o), 0)
})

test_that("a single population fails fast with a clear error", {
  # [approved F9] previously an opaque "subscript out of bounds"
  g1 <- gl.keep.pop(testset.gl, popNames(testset.gl)[1], verbose = 0)
  expect_error(gl.dist.pop(g1, method = "euclidean", plot.display = FALSE,
                           verbose = 0), "at least two")
})

test_that("SilicoDArT anchors and the sorensen->simple coercion", {
  gsp <- gsp4()
  anchors <- list(
    euclidean = c(5.051618, 3.957996, 4.90631, 4.060031, 5.52712,
                  4.94194),
    simple    = c(0.152543, 0.122633, 0.139212, 0.116141, 0.150073,
                  0.141082),
    jaccard   = c(0.395615, 0.322566, 0.37621, 0.306148, 0.396981,
                  0.367414)
  )
  for (m in names(anchors)) {
    D <- gl.dist.pop(gsp, method = m, plot.display = FALSE, verbose = 0)
    expect_equal(round(as.vector(D), 6), anchors[[m]], info = m)
  }
  # [approved: chain F2 -- completes when review-gl.dist.ind merges]
  # sorensen routes through gl.dist.ind. Before that branch merges,
  # gl.dist.ind coerces it to simple matching with a leaked warning;
  # after it merges, gl.dist.pop returns the real Sorensen collapse.
  # This assertion detects the chain state and pins the correct
  # behaviour for each.
  o <- capture.output(Dso <- gl.dist.pop(gsp, method = "sorensen",
                                         plot.display = FALSE,
                                         verbose = 0))
  if (length(o) > 0) {
    # pre-chain state (gl.dist.ind fix not yet merged)
    expect_equal(round(as.vector(Dso), 6), anchors$simple)
  } else {
    # post-chain state: real Sorensen through the collapse (via
    # gl.dist.ind, which carries the individual labels the collapse
    # requires -- the same route gl.dist.pop takes)
    Deng <- gl.dist.ind(gsp, method = "sorensen", plot.display = FALSE,
                        verbose = 0)
    Dsor <- utils.collapse.matrix(D = Deng, x = gsp, verbose = 0)
    expect_equal(as.vector(Dso), as.vector(Dsor))
    expect_false(isTRUE(all.equal(round(as.vector(Dso), 6),
                                  anchors$simple)))
  }
})

test_that("datatype gating: refused method/datatype pairs stop", {
  expect_error(gl.dist.pop(gsp4(), method = "nei", plot.display = FALSE,
                           verbose = 0))
  expect_error(gl.dist.pop(glp4(), method = "jaccard",
                           plot.display = FALSE, verbose = 0))
})

test_that("the input object is not modified", {
  glp <- glp4()
  h0 <- length(glp@other$history)
  invisible(gl.dist.pop(glp, method = "euclidean", plot.display = FALSE,
                        verbose = 0))
  expect_equal(length(glp@other$history), h0)
})
