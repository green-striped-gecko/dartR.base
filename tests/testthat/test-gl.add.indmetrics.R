# Characterization tests for gl.add.indmetrics
# Baseline snapshotted before review (review-gl.add.indmetrics).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect the approved behaviour change (the
# superset/partial-overlap metadata case now succeeds).

make_sub <- function(n = 6) {
  s <- testset.gl[1:n, ]
  s@other$loc.metrics <- testset.gl@other$loc.metrics
  s
}

test_that("gl.add.indmetrics: documented example (metadata matches the dart file)", {
  dartfile <- system.file("extdata", "testset_SNPs_2Row.csv",
                          package = "dartR.data")
  metadata <- system.file("extdata", "testset_metadata.csv",
                          package = "dartR.data")
  gl <- gl.read.dart(dartfile, probar = FALSE, verbose = 0)
  out <- capture.output(
    r <- gl.add.indmetrics(gl, ind.metafile = metadata, verbose = 0)
  )
  expect_equal(nInd(r), 250)
  expect_equal(nrow(r@other$ind.metrics), 250)
  expect_false(is.null(pop(r)))
})

test_that("gl.add.indmetrics keeps the matched individuals when metadata has extra rows", {
  # [approved diff F1] a metadata file whose ids are a superset of x now
  # succeeds (matched individuals returned with their metadata), where it
  # previously crashed at `ind.cov$pop_old <- x@pop`
  td <- tempdir()
  sub <- make_sub(6)
  ids <- indNames(sub)
  meta <- data.frame(
    id = c(ids, paste0("EXTRA_", 1:4)),
    pop = c(rep("P1", 6), rep("PX", 4)),
    lat = seq(-35, -33, length.out = 10),
    lon = seq(148, 150, length.out = 10),
    sex = rep(c("M", "F"), 5)
  )
  mf <- file.path(td, "im-extra.csv"); write.csv(meta, mf, row.names = FALSE)
  out <- capture.output(
    r <- gl.add.indmetrics(sub, ind.metafile = mf, verbose = 0)
  )
  expect_equal(nInd(r), 6)
  expect_equal(nrow(r@other$ind.metrics), 6)
  # only the 6 matched individuals, with their own metadata
  expect_setequal(indNames(r), ids)
  expect_true(all(as.character(pop(r)) == "P1"))
  expect_true(all(r@other$ind.metrics$sex %in% c("M", "F")))
})

test_that("gl.add.indmetrics aligns metadata to x for an equal-size reordered file", {
  td <- tempdir()
  sub <- make_sub(6)
  ids <- indNames(sub)
  meta <- data.frame(
    id = ids[c(3, 1, 5, 2, 6, 4)],
    pop = paste0("P", c(3, 1, 5, 2, 6, 4)),
    lat = seq(-35, -34, length.out = 6),
    lon = seq(148, 149, length.out = 6)
  )
  mf <- file.path(td, "im-reord.csv"); write.csv(meta, mf, row.names = FALSE)
  out <- capture.output(
    r <- gl.add.indmetrics(sub, ind.metafile = mf, verbose = 0)
  )
  expect_equal(nInd(r), 6)
  expect_equal(as.character(pop(r)),
               meta$pop[match(indNames(r), meta$id)])
  expect_true(is.numeric(r@other$latlon[, 1]))
})

test_that("gl.add.indmetrics subsets x to a strict metadata subset", {
  td <- tempdir()
  sub <- make_sub(6)
  ids <- indNames(sub)
  meta <- data.frame(id = ids[1:3], pop = c("A", "B", "C"),
                     lat = c(-34, -34.5, -35), lon = c(148, 148.5, 149))
  mf <- file.path(td, "im-subset.csv"); write.csv(meta, mf, row.names = FALSE)
  out <- capture.output(
    r <- gl.add.indmetrics(sub, ind.metafile = mf, verbose = 0)
  )
  expect_equal(nInd(r), 3)
  expect_equal(nrow(r@other$ind.metrics), 3)
})

test_that("gl.add.indmetrics is silent at verbose 0 and errors without an id column", {
  td <- tempdir()
  sub <- make_sub(6)
  ids <- indNames(sub)
  meta <- data.frame(id = ids, pop = "P", lat = -34, lon = 148)
  mf <- file.path(td, "im-ok.csv"); write.csv(meta, mf, row.names = FALSE)
  out <- capture.output(r <- gl.add.indmetrics(sub, ind.metafile = mf,
                                               verbose = 0))
  expect_length(out, 0)
  meta2 <- data.frame(name = ids, pop = "P")
  mf2 <- file.path(td, "im-noid.csv"); write.csv(meta2, mf2, row.names = FALSE)
  expect_error(gl.add.indmetrics(sub, ind.metafile = mf2, verbose = 0),
               "no id column")
})

test_that("gl.add.indmetrics errors clearly on duplicate ids", {
  # [approved diff F3] now stop(error(...)) with a message, not cat()+stop()
  td <- tempdir()
  sub <- make_sub(6)
  ids <- indNames(sub)
  meta <- data.frame(id = c(ids[1], ids[1], ids[3:6]), pop = "P")
  mf <- file.path(td, "im-dup.csv"); write.csv(meta, mf, row.names = FALSE)
  expect_error(gl.add.indmetrics(sub, ind.metafile = mf, verbose = 0),
               "not unique")
})

test_that("gl.add.indmetrics appends a history entry", {
  td <- tempdir()
  sub <- make_sub(6)
  meta <- data.frame(id = indNames(sub), pop = "P", lat = -34, lon = 148)
  mf <- file.path(td, "im-h.csv"); write.csv(meta, mf, row.names = FALSE)
  h0 <- length(sub@other$history)
  r <- gl.add.indmetrics(sub, ind.metafile = mf, verbose = 0)
  expect_equal(length(r@other$history), h0 + 1)
})
