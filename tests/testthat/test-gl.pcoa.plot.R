# Characterization tests for gl.pcoa.plot (function-review campaign).
# Written against upstream/dev ddaed27 to pin the pre-review behaviour, then
# updated for the changes approved on 2026-09-07. Every expectation that moved
# carries an "# [approved Fn]" comment naming the finding that moved it; every
# other expectation is unchanged from the baseline and must still hold.
# Finding IDs refer to function-review/reports/dartR.base/gl.pcoa.plot.md.
#
# Introspection is on the ggplot OBJECT (labels, layer data, mappings), not
# rendered images. Plotly branches are built with plotly_build(); nothing is
# sent to a viewer.

pdf(NULL)
op <- options(browser = function(...) invisible(NULL))
withr_restore <- function() options(op)

skip_if_not_installed("directlabels")

# Fixtures, built once per file run
.fx <- new.env()
fx <- function() {
  if (is.null(.fx$pca)) {
    .fx$gl <- dartR.data::testset.gl
    quiet <- capture.output({
      .fx$pca  <- gl.pcoa(.fx$gl, nfactors = 5, verbose = 0)
      .fx$pca2 <- gl.pcoa(.fx$gl, nfactors = 2, verbose = 0)
      .fx$gs30 <- .fx$gl[1:30, ]
      .fx$D    <- gl.dist.ind(.fx$gs30, verbose = 0)
      .fx$pco  <- gl.pcoa(.fx$D, nfactors = 3, verbose = 0)
    })
  }
  .fx
}
pct <- function(eig) round(eig * 100 / sum(eig[eig >= 0]), 1)

# ---------------------------------------------------------------- axis truth

test_that("2D axis labels state 100*eig/sum(positive eig) for the chosen axes", {
  f <- fx()
  e <- pct(f$pca$eig)
  quiet <- capture.output(p <- gl.pcoa.plot(f$pca, f$gl, verbose = 0))
  expect_equal(p$labels$x, paste("PCA Axis", 1, "(", e[1], "%)"))
  expect_equal(p$labels$y, paste("PCA Axis", 2, "(", e[2], "%)"))
  quiet <- capture.output(
    p25 <- gl.pcoa.plot(f$pca, f$gl, xaxis = 2, yaxis = 5, verbose = 0))
  expect_equal(p25$labels$x, paste("PCA Axis", 2, "(", e[2], "%)"))
  expect_equal(p25$labels$y, paste("PCA Axis", 5, "(", e[5], "%)"))
})

test_that("dist-derived ordination (no correction) is labelled PCoA with
           positive-eig denominator", {
  f <- fx()
  expect_null(f$pco$loadings)          # drives the PCoA classification
  expect_gt(sum(f$pco$eig < 0), 0)     # negative eigenvalues present
  e <- pct(f$pco$eig)
  quiet <- capture.output(pd <- gl.pcoa.plot(f$pco, f$gs30, verbose = 0))
  expect_equal(pd$labels$x, paste("PCoA Axis", 1, "(", e[1], "%)"))
  expect_equal(pd$labels$y, paste("PCoA Axis", 2, "(", e[2], "%)"))
})

test_that("dist-derived ordination WITH correction is labelled 'PCoA Axis'", {
  f <- fx()
  quiet <- capture.output(
    pcoc <- gl.pcoa(f$D, nfactors = 3, correction = "cailliez", verbose = 0))
  expect_false(is.null(pcoc$loadings)) # vectors.cor present in $loadings
  # Loadings are one row per entity, not one per locus -- the PCoA signal
  expect_equal(nrow(as.matrix(pcoc$loadings)), nrow(pcoc$scores))
  quiet <- capture.output(pdc <- gl.pcoa.plot(pcoc, f$gs30, verbose = 0))
  expect_match(pdc$labels$x, "^PCoA Axis 1")  # [approved F5] was "PCA Axis 1"
  e <- pct(pcoc$eig)
  expect_equal(pdc$labels$x, paste("PCoA Axis", 1, "(", e[1], "%)"))
  expect_equal(pdc$data$PCoAx, unname(pcoc$scores[, 1]))
})

# ---------------------------------------------------------------- data truth

test_that("plotted coordinates are the score columns; ind/pop map 1:1", {
  f <- fx()
  quiet <- capture.output(p <- gl.pcoa.plot(f$pca, f$gl, verbose = 0))
  expect_equal(p$data$PCoAx, unname(f$pca$scores[, 1]))
  expect_equal(p$data$PCoAy, unname(f$pca$scores[, 2]))
  expect_identical(as.character(p$data$ind), indNames(f$gl))
  expect_identical(as.character(p$data$pop), as.character(pop(f$gl)))
  quiet <- capture.output(
    p25 <- gl.pcoa.plot(f$pca, f$gl, xaxis = 2, yaxis = 5, verbose = 0))
  expect_equal(p25$data$PCoAx, unname(f$pca$scores[, 2]))
  expect_equal(p25$data$PCoAy, unname(f$pca$scores[, 5]))
})

test_that("a shuffled pop factor tracks the plot's pop mapping", {
  f <- fx()
  gl.sh <- f$gl
  set.seed(42)
  pop(gl.sh) <- sample(pop(f$gl))
  quiet <- capture.output(psh <- gl.pcoa.plot(f$pca, gl.sh, verbose = 0))
  expect_identical(as.character(psh$data$pop), as.character(pop(gl.sh)))
  expect_false(identical(as.character(psh$data$pop), as.character(pop(f$gl))))
})

test_that("as.pop substitutes the named ind.metric for pop", {
  f <- fx()
  quiet <- capture.output(
    pap <- gl.pcoa.plot(f$pca, f$gl, as.pop = "sex", verbose = 0))
  expect_identical(as.character(pap$data$pop),
                   as.character(f$gl@other$ind.metrics$sex))
  quiet <- capture.output(
    err <- tryCatch(gl.pcoa.plot(f$pca, f$gl, as.pop = "nonexistent",
                                 verbose = 0),
                    error = function(e) conditionMessage(e)))
  expect_match(err, "ind.metrics")  # [approved F14] was "loc.metrics"
})

# ------------------------------------------------------- layers and aesthetics

test_that("layer census: pop branch 4 layers, none branch 3, ellipse appends
           stat_ellipse at plevel", {
  f <- fx()
  quiet <- capture.output(p <- gl.pcoa.plot(f$pca, f$gl, verbose = 0))
  expect_length(p$layers, 4)
  expect_identical(unname(sapply(p$layers, function(l) class(l$geom)[1])),
                   c("GeomPoint", "GeomDl", "GeomHline", "GeomVline"))
  quiet <- capture.output(
    pn <- gl.pcoa.plot(f$pca, f$gl, pop.labels = "none", verbose = 0))
  expect_length(pn$layers, 3)
  quiet <- capture.output(
    pe <- gl.pcoa.plot(f$pca, f$gl, ellipse = TRUE, plevel = 0.9, verbose = 0))
  expect_length(pe$layers, 5)
  lev <- pe$layers[[5]]$stat_params$level
  expect_equal(lev, 0.9)
  expect_equal(pe$layers[[5]]$stat_params$type, "norm")
})

test_that("pt.size, pt.colors, pt.shapes reach the plot object (2D pop branch)", {
  f <- fx()
  quiet <- capture.output(pz <- gl.pcoa.plot(f$pca, f$gl, pt.size = 5, verbose = 0))
  expect_equal(pz$layers[[1]]$aes_params$size, 5)
  npop <- nlevels(pop(f$gl))
  cols <- grDevices::rainbow(npop)
  shps <- rep(c(16, 17, 15, 0, 2), length.out = npop)
  quiet <- capture.output(
    pc2 <- gl.pcoa.plot(f$pca, f$gl, pt.colors = cols, pt.shapes = shps,
                        verbose = 0))
  bld <- ggplot2::ggplot_build(pc2)
  used <- cols[as.integer(factor(pop(f$gl)))]
  expect_setequal(unique(bld$data[[1]]$colour), unique(used))
  expect_setequal(unique(bld$data[[1]]$shape),
                  unique(shps[as.integer(factor(pop(f$gl)))]))
})

test_that("hadjust and vadjust set the label justification [approved F2]", {
  f <- fx()
  quiet <- capture.output(
    p0 <- gl.pcoa.plot(f$pca, f$gl, hadjust = 0, vadjust = 0, verbose = 0))
  quiet <- capture.output(
    p3 <- gl.pcoa.plot(f$pca, f$gl, hadjust = 3, vadjust = 3, verbose = 0))
  d0 <- ggplot2::ggplot_build(p0)$data[[2]]   # the GeomDl label layer
  d3 <- ggplot2::ggplot_build(p3)$data[[2]]
  # [approved F2] the baseline asserted these two builds were identical
  expect_false(identical(d0, d3))
  expect_true(all(d0$hjust == 0) && all(d0$vjust == 0))
  expect_true(all(d3$hjust == 3) && all(d3$vjust == 3))
  # the points themselves are untouched by the label justification
  expect_equal(ggplot2::ggplot_build(p0)$data[[1]],
               ggplot2::ggplot_build(p3)$data[[1]])
  # the range check now tests vadjust, and resets it to the documented 1
  out <- capture.output(
    pv <- gl.pcoa.plot(f$pca, f$gl, vadjust = 5, verbose = 2))
  expect_gt(length(grep("vadjust", out)), 0)   # [approved F10] now warns
  expect_length(grep("hadjust", out), 0)
  expect_true(all(ggplot2::ggplot_build(pv)$data[[2]]$vjust == 1))
})

test_that("scale = TRUE fixes the coordinate ratio at exactly 1 for any axis
           pair [F15]", {
  f <- fx()
  quiet <- capture.output(
    ps <- gl.pcoa.plot(f$pca, f$gl, scale = TRUE, verbose = 0))
  expect_equal(ggplot2::ggplot_build(ps)$layout$coord$ratio, 1)
  quiet <- capture.output(
    ps15 <- gl.pcoa.plot(f$pca, f$gl, scale = TRUE, xaxis = 1, yaxis = 5,
                         verbose = 0))
  expect_equal(ggplot2::ggplot_build(ps15)$layout$coord$ratio, 1)  # [F15] not e-based
  quiet <- capture.output(
    pns <- gl.pcoa.plot(f$pca, f$gl, scale = FALSE, verbose = 0))
  expect_null(ggplot2::ggplot_build(pns)$layout$coord$ratio)
})

test_that("legend branch titles the colour legend 'Population' [approved F16]", {
  f <- fx()
  quiet <- capture.output(
    pg <- gl.pcoa.plot(f$pca, f$gl, pop.labels = "legend", verbose = 0))
  bg <- ggplot2::ggplot_build(pg)
  expect_equal(as.character(bg$plot$labels$colour), "Population")  # was "pop"
  # the colours themselves are unchanged by dropping the layer-level aes
  expect_identical(as.character(pg$data$pop), as.character(pop(f$gl)))
})

# ------------------------------------------------------------- branch contract

test_that("returns the plot object visibly, as @return now documents [F8]", {
  f <- fx()
  quiet <- capture.output(
    res <- withVisible(gl.pcoa.plot(f$pca, f$gl, verbose = 0)))
  expect_s3_class(res$value, "ggplot")
  expect_true(res$visible)  # unchanged; @return was corrected to match
})

test_that("interactive branch returns a plotly object, silent at verbose = 0
           [approved F6, F13]", {
  skip_if_not_installed("plotly")
  f <- fx()
  out <- capture.output(
    pi2 <- gl.pcoa.plot(f$pca, f$gl, interactive = TRUE, verbose = 0))
  expect_s3_class(pi2, "plotly")
  expect_length(out, 0)             # [approved F6] baseline leaked 3 lines
  # the false NOTE is gone; what is printed at verbose 2 is truthful [F13]
  out2 <- capture.output(
    pi3 <- gl.pcoa.plot(f$pca, f$gl, interactive = TRUE, verbose = 2))
  expect_length(grep("Returning the ordination scores", out2), 0)
  expect_gt(length(grep("Returning an interactive plotly object", out2)), 0)
})

test_that("interactive branch honours pt.colors and pt.shapes [approved F12]", {
  skip_if_not_installed("plotly")
  f <- fx()
  npop <- nlevels(pop(f$gl))
  cols <- grDevices::rainbow(npop)
  shps <- rep(c(16, 17, 15, 0, 2), length.out = npop)
  quiet <- capture.output(
    pin <- gl.pcoa.plot(f$pca, f$gl, interactive = TRUE, pt.colors = cols,
                        pt.shapes = shps, verbose = 0))
  bb <- suppressWarnings(plotly::plotly_build(pin))
  rgba2hex <- function(s) {
    v <- as.numeric(strsplit(gsub("rgba?\\(|\\)", "", s), ",")[[1]])
    toupper(grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255))
  }
  used <- unique(unlist(lapply(bb$x$data, function(tr) tr$marker$color)))
  used <- unique(vapply(used, rgba2hex, ""))
  expect_gt(length(used), 1)
  expect_true(all(used %in% toupper(cols)))   # baseline: plotly default palette
  sym <- unique(unlist(lapply(bb$x$data, function(tr) tr$marker$symbol)))
  expect_setequal(sym, c("circle", "triangle-up", "square",
                         "square-open", "triangle-up-open"))
})

test_that("3D branch returns plotly with truthful axis labels and score
           coordinates", {
  skip_if_not_installed("plotly")
  f <- fx()
  e <- pct(f$pca$eig)
  out <- capture.output(
    p3 <- gl.pcoa.plot(f$pca, f$gl, zaxis = 3, verbose = 0))
  expect_s3_class(p3, "plotly")
  expect_length(out, 0)
  bb <- suppressWarnings(plotly::plotly_build(p3))
  sc <- bb$x$layout$scene
  expect_equal(sc$xaxis$title, paste("PCA Axis", 1, "(", e[1], "%)"))
  expect_equal(sc$yaxis$title, paste("PCA Axis", 2, "(", e[2], "%)"))
  expect_equal(sc$zaxis$title, paste("PCA Axis", 3, "(", e[3], "%)"))
  xs <- sort(unlist(lapply(bb$x$data, function(tr) tr$x)))
  expect_equal(xs, sort(unname(f$pca$scores[, 1])))
})

test_that("SilicoDArT ordination plots with truthful labels and coordinates", {
  f <- fx()
  quiet <- capture.output({
    pgs  <- gl.pcoa(dartR.data::testset.gs, nfactors = 3, verbose = 0)
    psil <- gl.pcoa.plot(pgs, dartR.data::testset.gs, verbose = 0)
  })
  e <- pct(pgs$eig)
  expect_equal(psil$labels$x, paste("PCA Axis", 1, "(", e[1], "%)"))
  expect_equal(psil$data$PCoAx, unname(pgs$scores[, 1]))
})

# ------------------------------------------------------------------ error paths

test_that("pop.labels = 'ind' builds a plot labelled with individual names
           [approved F1]", {
  f <- fx()
  quiet <- capture.output(
    pind <- gl.pcoa.plot(f$pca, f$gl, pop.labels = "ind", verbose = 0))
  # [approved F1] the baseline crashed here with "object 'plott' not found"
  expect_s3_class(pind, "ggplot")
  expect_length(pind$layers, 4)
  expect_identical(unname(sapply(pind$layers, function(l) class(l$geom)[1])),
                   c("GeomPoint", "GeomDl", "GeomHline", "GeomVline"))
  b <- ggplot2::ggplot_build(pind)
  expect_identical(as.character(b$data[[2]]$label), indNames(f$gl))
  # data truth holds on the new branch
  expect_equal(pind$data$PCoAx, unname(f$pca$scores[, 1]))
  expect_equal(pind$data$PCoAy, unname(f$pca$scores[, 2]))
})

test_that("axis bounds are checked, not overwritten with constants
           [approved F3]", {
  f <- fx()
  quiet <- capture.output(
    pca1 <- gl.pcoa(f$gl, nfactors = 1, verbose = 0))
  quiet <- capture.output(
    err1 <- tryCatch(gl.pcoa.plot(pca1, f$gl, verbose = 0),
                     error = function(e) conditionMessage(e)))
  # [approved F3] baseline: "subscript out of bounds"
  expect_no_match(err1, "subscript out of bounds")
  expect_match(err1, "at least 2 axes")
  quiet <- capture.output(
    err2 <- tryCatch(gl.pcoa.plot(f$pca2, f$gl, zaxis = 5, verbose = 0),
                     error = function(e) conditionMessage(e)))
  expect_no_match(err2, "subscript out of bounds")
  expect_match(err2, "at least 3 axes")
  # a two-individual object yields a one-axis ordination: same clear message
  gl2 <- f$gl[1:2, ]
  quiet <- capture.output(pp2 <- gl.pcoa(gl2, nfactors = 2, verbose = 0))
  quiet <- capture.output(
    err3 <- tryCatch(gl.pcoa.plot(pp2, gl2, verbose = 0),
                     error = function(e) conditionMessage(e)))
  expect_no_match(err3, "subscript out of bounds")
  expect_match(err3, "at least 2 axes")
  # the chosen axes must differ from one another
  quiet <- capture.output(
    err4 <- tryCatch(gl.pcoa.plot(f$pca, f$gl, xaxis = 2, yaxis = 2,
                                  verbose = 0),
                     error = function(e) conditionMessage(e)))
  expect_match(err4, "must differ")
})

test_that("in-range axis clamp warnings are gated at verbose 2 [approved F6]", {
  f <- fx()
  out <- capture.output(
    zx <- gl.pcoa.plot(f$pca, f$gl, xaxis = 9, verbose = 0))
  expect_length(grep("X-axis", out), 0)      # [approved F6] baseline printed
  expect_match(zx$labels$x, "^PCA Axis 1")   # still clamped to 1
  out.v2 <- capture.output(
    zx2 <- gl.pcoa.plot(f$pca, f$gl, xaxis = 9, verbose = 2))
  expect_gt(length(grep("X-axis", out.v2)), 0)
  out2 <- capture.output(
    zp <- gl.pcoa.plot(f$pca, f$gl, plevel = 2, verbose = 0))
  expect_length(grep("plevel", out2), 0)     # [approved F6] baseline printed
})

test_that("unpaired inputs fail with a message that names the problem
           [approved F11]", {
  f <- fx()
  quiet <- capture.output({
    Dp   <- gl.dist.pop(f$gl, verbose = 0)
    pcop <- gl.pcoa(Dp, nfactors = 3, verbose = 0)
  })
  quiet <- capture.output(
    err <- tryCatch(gl.pcoa.plot(pcop, f$gl, verbose = 0),
                    error = function(e) conditionMessage(e)))
  # [approved F11] baseline: "arguments imply differing number of rows: 30, 250"
  expect_no_match(err, "differing number of rows")
  expect_match(err, "entities but the genlight object holds")
  expect_match(err, "gl.dist.pop")
  # a data.frame classifies as "list"; it is now rejected before the
  # animation branch can subset a genlight [approved F11]
  quiet <- capture.output(
    err2 <- tryCatch(gl.pcoa.plot(as.data.frame(f$pca$scores), f$gl,
                                  verbose = 0),
                     error = function(e) conditionMessage(e)))
  expect_no_match(err2, "not subsettable")
  expect_match(err2, "animation plot requires")
})

test_that("plainly wrong classes are rejected informatively by the gate", {
  f <- fx()
  quiet <- capture.output(
    err <- tryCatch(gl.pcoa.plot(f$pca$scores, f$gl, verbose = 0),
                    error = function(e) conditionMessage(e)))
  expect_match(err, "found matrix expecting glPca or list")
  quiet <- capture.output(
    err2 <- tryCatch(gl.pcoa.plot(f$pco, f$D, verbose = 0),
                     error = function(e) conditionMessage(e)))
  expect_match(err2, "found dist expecting SNP or SilicoDArT or fd or list")
})

# ------------------------------------------------------------ side effects, IO

test_that("verbose = 0 is silent in text and draws nothing [approved F6, F7]", {
  f <- fx()
  out <- capture.output(p <- gl.pcoa.plot(f$pca, f$gl, verbose = 0))
  expect_length(out, 0)
  out2 <- capture.output(
    pn <- gl.pcoa.plot(f$pca, f$gl, pop.labels = "none", verbose = 0))
  expect_length(out2, 0)  # [approved F6] baseline printed the "none" message
  # VRB5: silence covers graphics too. The device display list stays empty
  # unless something was drawn.
  drew <- function(expr) {
    fx2 <- tempfile(fileext = ".pdf")
    grDevices::pdf(fx2)
    grDevices::dev.control(displaylist = "enable")
    on.exit({grDevices::dev.off(); unlink(fx2)}, add = TRUE)
    force(expr)
    length(grDevices::recordPlot()[[1]]) > 0
  }
  quiet <- capture.output(d0 <- drew(gl.pcoa.plot(f$pca, f$gl, verbose = 0)))
  expect_false(d0)                      # [approved F7] baseline always drew
  quiet <- capture.output(d2 <- drew(gl.pcoa.plot(f$pca, f$gl, verbose = 2)))
  expect_true(d2)                       # plot.display defaults to TRUE
  quiet <- capture.output(
    d3 <- drew(gl.pcoa.plot(f$pca, f$gl, plot.display = FALSE, verbose = 2)))
  expect_false(d3)
  # PLT3: the returned object does not depend on the display
  quiet <- capture.output(
    pd <- gl.pcoa.plot(f$pca, f$gl, plot.display = FALSE, verbose = 0))
  expect_s3_class(pd, "ggplot")
  expect_equal(pd$data$PCoAx, unname(f$pca$scores[, 1]))
})

test_that("inputs come back untouched; no cwd or sink drift", {
  f <- fx()
  gl.copy <- f$gl
  pca.copy <- f$pca
  wd <- getwd(); sn <- sink.number()
  quiet <- capture.output(tmp <- gl.pcoa.plot(f$pca, gl.copy, verbose = 0))
  expect_identical(gl.copy, f$gl)
  expect_identical(pca.copy, f$pca)
  expect_identical(getwd(), wd)
  expect_identical(sink.number(), sn)
})

test_that("plot.file saves an RDS to plot.dir; with plot.dir NULL it lands in
           tempdir, not getwd() [approved F4]", {
  f <- fx()
  dir1 <- file.path(tempdir(), "pcoaplot-explicit")
  dir.create(dir1, showWarnings = FALSE)
  quiet <- capture.output(
    z <- gl.pcoa.plot(f$pca, f$gl, plot.file = "chk", plot.dir = dir1,
                      verbose = 0))
  expect_true(file.exists(file.path(dir1, "chk.RDS")))
  expect_s3_class(readRDS(file.path(dir1, "chk.RDS")), "ggplot")
  dir2 <- file.path(tempdir(), "pcoaplot-cwd")
  dir.create(dir2, showWarnings = FALSE)
  old <- setwd(dir2)
  on.exit(setwd(old), add = TRUE)
  quiet <- capture.output(
    z <- gl.pcoa.plot(f$pca, f$gl, plot.file = "cwdplot", verbose = 0))
  setwd(old)
  # [approved F4] baseline wrote cwdplot.RDS into the working directory
  expect_false(file.exists(file.path(dir2, "cwdplot.RDS")))
  expect_true(file.exists(file.path(tempdir(), "cwdplot.RDS")))
})

withr_restore()
