# Characterization tests for gl.map.interactive
# Baseline snapshotted before review (review-gl.map.interactive) at be7cb2e.
# Assertions marked [approved diff] were flipped in Phase C to reflect the
# approved behaviour changes (report: function-review/reports/dartR.base/
# gl.map.interactive.md). Introspection is on the leaflet widget's call list
# (m$x$calls), not on rendered output.

skip_if_not_installed("leaflet")
skip_if_not_installed("leaflet.minicharts")

map_calls <- function(m) vapply(m$x$calls, function(z) z$method, character(1))
call_args <- function(m, method, which = 1) {
  m$x$calls[[which(map_calls(m) == method)[which]]]$args
}
# leaflet encodes addCircles as (lat, lng, radius, layerId, group, options,
# popup, ...) and addMarkers as (lat, lng, icon, layerId, group, options,
# popup, popupOptions, clusterOptions, clusterId, label, labelOptions, ...)
polyline_segments <- function(m) {
  lapply(m$x$calls[map_calls(m) == "addPolylines"], function(p) {
    seg <- p$args[[1]][[1]][[1]][[1]]
    list(lng = seg$lng, lat = seg$lat, color = p$args[[4]]$color,
         weight = p$args[[4]]$weight)
  })
}
flow_colors <- function(m) {
  vapply(m$x$calls[map_calls(m) == "addFlows"],
         function(p) p$args[[1]][[1]]$static$color, character(1))
}
seg_ends <- function(s, ll) sort(match(round(s$lng, 6), round(ll$lon, 6)))

small_platypus <- function(n = 6) {
  x <- platypus.gl[1:n, ]
  x@other$latlon <- platypus.gl@other$latlon[1:n, ]
  x@other$ind.metrics <- platypus.gl@other$ind.metrics[1:n, ]
  x
}

test_that("default call on platypus.gl: widget structure and data truth", {
  x <- platypus.gl
  out <- capture.output(m <- gl.map.interactive(x, verbose = 0))
  expect_length(out, 0)
  expect_s3_class(m, "leaflet")
  # [approved diff, change 10] one base layer instead of addTiles + provider
  expect_equal(map_calls(m), c("addProviderTiles", "addCircles", "addMarkers",
                               "addScaleBar"))
  ci <- call_args(m, "addCircles")
  expect_equal(unlist(ci[[1]]), x@other$latlon$lat)
  expect_equal(unlist(ci[[2]]), x@other$latlon$lon)
  expect_equal(ci[[6]]$weight, 10)
  expect_equal(ci[[6]]$opacity, 0.8)
  expect_equal(length(unique(ci[[6]]$color)), nPop(x))
  expect_equal(unlist(ci[[7]]), indNames(x))
  mk <- call_args(m, "addMarkers")
  expect_equal(unlist(mk[[11]]), popNames(x))
  expect_equal(mk[[12]]$textsize, "12px")
  cen <- apply(x@other$latlon, 2, function(v) tapply(v, pop(x), mean))
  expect_equal(unname(unlist(mk[[1]])), unname(cen[, "lat"]))
  expect_equal(unname(unlist(mk[[2]])), unname(cen[, "lon"]))
  expect_equal(call_args(m, "addProviderTiles")[[1]], "Esri.NatGeoWorldMap")
  # read-only: input untouched
  expect_identical(x, platypus.gl)
})

test_that("switches: no circles, no labels, no scale bar, provider", {
  m <- gl.map.interactive(platypus.gl, ind.circles = FALSE, pop.labels = FALSE,
                          scale.bar = FALSE, provider = "OpenStreetMap",
                          verbose = 0)
  expect_equal(map_calls(m), "addProviderTiles")   # [approved diff, change 10]
  expect_equal(call_args(m, "addProviderTiles")[[1]], "OpenStreetMap")
})

test_that("SilicoDArT and FBM objects are accepted", {
  m <- gl.map.interactive(testset.gs, verbose = 0)
  expect_equal(length(unlist(call_args(m, "addCircles")[[1]])), nInd(testset.gs))
  xf <- gl.gen2fbm(platypus.gl, verbose = 0)
  m <- gl.map.interactive(xf, verbose = 0)
  expect_s3_class(m, "leaflet")
})

test_that("F1: single population label sits at the centre for either column order", {
  x1 <- gl.keep.pop(platypus.gl, pop.list = popNames(platypus.gl)[1],
                    verbose = 0)
  m <- gl.map.interactive(x1, verbose = 0)
  mk <- call_args(m, "addMarkers")
  # [approved diff, change 1] was lat <- mean(lon), lng <- mean(lat)
  expect_equal(unname(unlist(mk[[1]])), mean(x1@other$latlon$lat))
  expect_equal(unname(unlist(mk[[2]])), mean(x1@other$latlon$lon))
  xb <- gl.keep.pop(bandicoot.gl, pop.list = popNames(bandicoot.gl)[1],
                    verbose = 0)
  mb <- gl.map.interactive(xb, verbose = 0)
  mk <- call_args(mb, "addMarkers")
  expect_equal(unname(unlist(mk[[1]])), mean(xb@other$latlon$lat))
  expect_equal(unname(unlist(mk[[2]])), mean(xb@other$latlon$lon))
})

test_that("F2: dist objects from gl.dist.pop / gl.dist.ind are accepted", {
  d <- gl.dist.pop(platypus.gl, verbose = 0)
  expect_s3_class(d, "dist")
  # [approved diff, change 2] was "argument is of length zero"
  m <- gl.map.interactive(platypus.gl, matrix = d, verbose = 0)
  expect_equal(sum(map_calls(m) == "addPolylines"), 3)
  x <- small_platypus()
  di <- gl.dist.ind(x, verbose = 0)
  m <- gl.map.interactive(x, matrix = di, verbose = 0)
  expect_equal(sum(map_calls(m) == "addPolylines"), 15)
  expect_error(gl.map.interactive(x, matrix = "a", verbose = 0),
               "square matrix")
  expect_error(gl.map.interactive(x, matrix = matrix(0, 6, 3), verbose = 0),
               "square matrix")
})

test_that("population matrix, symmetric: one polyline per pair, none to self", {
  x <- platypus.gl
  dm <- as.matrix(gl.dist.pop(x, verbose = 0))
  out <- capture.output(m <- gl.map.interactive(x, matrix = dm, verbose = 0))
  expect_length(out, 0)   # [approved diff, change 5] gl.colors silenced
  segs <- polyline_segments(m)
  # [approved diff, change 4] was 6 lines, 3 of them self links
  expect_length(segs, 3)
  expect_equal(sum(vapply(segs, function(s) s$lng[1] == s$lng[2], logical(1))), 0)
  expect_true("addLegend" %in% map_calls(m))
  cen <- apply(x@other$latlon, 2, function(v) tapply(v, pop(x), mean))
  expect_equal(sort(unique(unlist(lapply(segs, `[[`, "lng")))),
               sort(unname(cen[, "lon"])))
  # [approved diff, change 9] line width follows the standardised value
  w <- vapply(segs, `[[`, numeric(1), "weight")
  expect_equal(sort(w), sort(unique(as.vector(
    (dm - min(dm)) / (max(dm) - min(dm)) * 9 + 1))[-1]), tolerance = 1e-8)
  expect_equal(range(w), c(1, 10) * 0 + range(w))
  lg <- call_args(m, "addLegend")[[1]]
  expect_equal(as.character(lg$labels), as.character(1:10))
  expect_equal(lg$group, "addPolylines")
  # standard = FALSE: raw values as width, zero pairs (diagonal) still skipped
  m0 <- gl.map.interactive(x, matrix = dm, standard = FALSE, verbose = 0)
  w0 <- vapply(polyline_segments(m0), `[[`, numeric(1), "weight")
  expect_equal(sort(w0), sort(dm[lower.tri(dm)]), tolerance = 1e-8)
})

test_that("F3: matrices are aligned by dimnames, else taken in object order", {
  x <- small_platypus()
  expect_equal(order(indNames(x)), c(4, 3, 5, 6, 1, 2))
  mm <- matrix(0, 6, 6, dimnames = list(indNames(x), indNames(x)))
  mm[2, 1] <- mm[1, 2] <- 5
  m <- gl.map.interactive(x, matrix = mm, standard = FALSE, verbose = 0)
  segs <- polyline_segments(m)
  ll <- x@other$latlon
  # [approved diff, changes 3 and 4] was 21 lines, hot ones at 1-6 and 2-5
  expect_length(segs, 1)
  expect_equal(seg_ends(segs[[1]], ll), c(1, 2))
  # rows/cols permuted but named: same link
  perm <- c(3, 1, 6, 2, 5, 4)
  mp <- mm[perm, perm]
  mperm <- gl.map.interactive(x, matrix = mp, standard = FALSE, verbose = 0)
  expect_equal(seg_ends(polyline_segments(mperm)[[1]], ll), c(1, 2))
  # unnamed: object order, no reorder
  mu <- unname(mm)
  munn <- gl.map.interactive(x, matrix = mu, standard = FALSE, verbose = 0)
  expect_equal(seg_ends(polyline_segments(munn)[[1]], ll), c(1, 2))
  # named but not matching: warning at verbose 2, object order
  mx <- mm; dimnames(mx) <- list(letters[1:6], letters[1:6])
  out <- capture.output(mnm <- gl.map.interactive(x, matrix = mx,
                                                  standard = FALSE, verbose = 2))
  expect_true(any(grepl("dimnames do not match", out)))
  expect_equal(seg_ends(polyline_segments(mnm)[[1]], ll), c(1, 2))
  # population matrix in reversed popNames order is realigned
  dm <- as.matrix(gl.dist.pop(platypus.gl, verbose = 0))
  ref <- polyline_segments(gl.map.interactive(platypus.gl, matrix = dm,
                                              verbose = 0))
  rev <- polyline_segments(gl.map.interactive(platypus.gl, matrix = dm[3:1, 3:1],
                                              verbose = 0))
  key <- function(s) paste(c(sort(round(s$lng, 6)), round(s$weight, 6)),
                           collapse = "_")
  expect_setequal(vapply(ref, key, ""), vapply(rev, key, ""))
})

test_that("asymmetric matrix draws flows; F6: NA cells are skipped or greyed", {
  x <- small_platypus()
  ma <- matrix(0, 6, 6, dimnames = list(indNames(x), indNames(x)))
  ma[2, 1] <- 5; ma[1, 2] <- 2
  m <- gl.map.interactive(x, matrix = ma, symmetric = FALSE,
                          standard = FALSE, verbose = 0)
  expect_equal(sum(map_calls(m) == "addFlows"), 30)
  expect_false("addLegend" %in% map_calls(m))
  expect_equal(sort(unique(flow_colors(m))), c("#00AA00", "#00AAFF", "#FFAA00"))
  ma[3, 4] <- NA
  # [approved diff, change 6] was "missing value where TRUE/FALSE needed"
  mn <- gl.map.interactive(x, matrix = ma, symmetric = FALSE,
                           standard = FALSE, verbose = 0)
  expect_equal(sum(map_calls(mn) == "addFlows"), 29)
  expect_equal(sum(flow_colors(mn) == "#333333"), 1)
})

test_that("F7: colour vector shorter than nPop errors", {
  # [approved diff, change 7] was NA colours for 76 of 81 individuals
  expect_error(gl.map.interactive(platypus.gl, ind.circle.cols = "red",
                                  verbose = 0), "colours but the dataset has")
  m <- gl.map.interactive(platypus.gl, ind.circle.cols = c("red", "blue", "green"),
                          verbose = 0)
  expect_setequal(unique(call_args(m, "addCircles")[[6]]$color),
                  c("red", "blue", "green"))
})

test_that("F8: latlon stored as a matrix is accepted", {
  xm <- platypus.gl
  xm@other$latlon <- as.matrix(platypus.gl@other$latlon)
  # [approved diff, change 7] was "$ operator is invalid for atomic vectors"
  m <- gl.map.interactive(xm, verbose = 0)
  expect_equal(unlist(call_args(m, "addCircles")[[2]]), platypus.gl@other$latlon$lon)
})

test_that("input validation", {
  xw <- platypus.gl
  colnames(xw@other$latlon) <- c("latitude", "longitude")
  expect_error(gl.map.interactive(xw, verbose = 0), "not named")
  xn <- platypus.gl
  xn@other$latlon <- NULL
  expect_error(gl.map.interactive(xn, verbose = 0), "No valid coordinates")
  expect_error(gl.map.interactive(platypus.gl, matrix = matrix(0, 4, 4),
                                  verbose = 0), "does neither match")
})
