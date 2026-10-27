# Per-genotype matrices in @other$geno.metrics (individuals x loci, e.g. the
# read counts and likelihoods kept by gl.read.vcf(read.counts = TRUE,
# likelihoods = TRUE)) must stay aligned with the genotypes whenever
# individuals or loci are removed, reordered or combined.

# Each cell codes its own individual and locus, so any misalignment changes
# the values, not only the dimensions
add_gm <- function(x) {
  ref <- outer(seq_len(nInd(x)) * 10000L, seq_len(nLoc(x)), "+")
  storage.mode(ref) <- "integer"
  x@other$geno.metrics <- list(ref.count = ref, alt.count = -ref,
                               pl.ref = ref + 1000000L,
                               pl.het = ref + 2000000L,
                               pl.alt = ref + 3000000L)
  x
}

# The matrices are kept by position, so the expected values are looked up
# by name in a labelled copy of the full matrices; ind gives the source
# individual of each row when the function renamed them
expect_aligned <- function(y, full, ind = indNames(y)) {
  gm <- y@other$geno.metrics
  expect_false(is.null(gm))
  for (nm in names(full@other$geno.metrics)) {
    m <- full@other$geno.metrics[[nm]]
    dimnames(m) <- list(indNames(full), locNames(full))
    expect_identical(gm[[nm]], unname(m[ind, locNames(y), drop = FALSE]))
  }
}

make_gm <- function() {
  add_gm(gl.compliance.check(testset.gl[1:20, 1:40], verbose = 0))
}

test_that("[ subsets geno.metrics by individuals and loci", {
  x <- make_gm()
  expect_aligned(x[c(3, 1, 7), c(5, 2, 9)], x)
  expect_aligned(x[, 10:1], x)
  expect_aligned(x[-(1:4), ], x)
  expect_aligned(x[c(TRUE, FALSE), locNames(x)[c(4, 8)]], x)
  # loci are dropped from the matrices even when other is not treated
  expect_aligned(x[1:5, 1:6, treatOther = FALSE], x)
})

test_that("[ keeps geno.metrics aligned when nInd equals 2 or nLoc", {
  x <- make_gm()
  # a list of two matrices has length 2: with two individuals the
  # per-individual pass would subset the list itself
  x2 <- x[1:2, ]
  expect_aligned(x2[2, 3:5], x)
  # with nInd == nLoc a dimension check alone cannot tell rows from columns
  sq <- x[, 1:20]
  expect_aligned(sq[c(2, 4), 5:1], x)
})

test_that("[ subsets geno.metrics on FBM-backed objects", {
  x <- gl.gen2fbm(make_gm(), verbose = 0)
  expect_aligned(x[c(3, 1), c(5, 2, 9)], x)
  expect_aligned(x[-1, -2], x)
})

test_that("rbind stacks geno.metrics in the first object's locus order", {
  x <- make_gm()
  r <- rbind(x[1:8, ], x[9:20, ])
  expect_aligned(r, x)
  r2 <- rbind(x[1:8, ], x[9:20, 40:1])
  expect_identical(locNames(r2), locNames(x))
  expect_aligned(r2, x)
})

test_that("cbind binds geno.metrics by loci", {
  x <- make_gm()
  expect_aligned(cbind(x[, 1:15], x[, 16:40]), x)
})

test_that("binding drops geno.metrics with a warning when an object lacks them", {
  x <- make_gm()
  y <- x[9:20, ]
  y@other$geno.metrics <- NULL
  expect_warning(r <- rbind(x[1:8, ], y), "geno.metrics")
  expect_null(r@other$geno.metrics)
  z <- x[, 16:40]
  z@other$geno.metrics <- NULL
  expect_warning(cc <- cbind(x[, 1:15], z), "geno.metrics")
  expect_null(cc@other$geno.metrics)
})

test_that("gl.compliance.check drops geno.metrics that no longer fit the object", {
  x <- make_gm()
  x@other$geno.metrics$ref.count <- x@other$geno.metrics$ref.count[, 1:39]
  out <- capture.output(y <- gl.compliance.check(x, verbose = 1))
  expect_null(y@other$geno.metrics)
  expect_true(any(grepl("geno.metrics", out)))
  # matrices that fit are left untouched
  x <- make_gm()
  capture.output(y <- gl.compliance.check(x, verbose = 0))
  expect_identical(y@other$geno.metrics, x@other$geno.metrics)
})

test_that("gl.randomize.snps swaps the read counts and likelihoods of the loci it recodes", {
  x <- make_gm()
  set.seed(7)
  y <- gl.randomize.snps(x, plot.display = FALSE, verbose = 0)
  swapped <- x@loc.all != y@loc.all
  expect_true(any(swapped) && !all(swapped))
  gx <- x@other$geno.metrics
  gy <- y@other$geno.metrics
  expect_identical(gy$ref.count[, swapped], gx$alt.count[, swapped])
  expect_identical(gy$alt.count[, swapped], gx$ref.count[, swapped])
  expect_identical(gy$ref.count[, !swapped], gx$ref.count[, !swapped])
  # the homozygote likelihoods swap; the heterozygote's stays
  expect_identical(gy$pl.ref[, swapped], gx$pl.alt[, swapped])
  expect_identical(gy$pl.alt[, swapped], gx$pl.ref[, swapped])
  expect_identical(gy$pl.ref[, !swapped], gx$pl.ref[, !swapped])
  expect_identical(gy$pl.het, gx$pl.het)
})

test_that("keep, drop, filter and subsample functions keep geno.metrics aligned", {
  x <- add_gm(gl.compliance.check(testset.gl[1:60, 1:200], verbose = 0))
  pops <- levels(pop(x))
  quiet <- function(expr) {
    capture.output(res <- suppressWarnings(expr))
    res
  }
  ys <- list(
    keep.ind = quiet(gl.keep.ind(x, indNames(x)[c(5, 2, 9)], verbose = 0)),
    drop.ind = quiet(gl.drop.ind(x, indNames(x)[1:3], verbose = 0)),
    keep.loc = quiet(gl.keep.loc(x, locNames(x)[c(30, 4, 7)], verbose = 0)),
    drop.loc = quiet(gl.drop.loc(x, locNames(x)[1:10], verbose = 0)),
    keep.pop = quiet(gl.keep.pop(x, pops[1], verbose = 0)),
    drop.pop = quiet(gl.drop.pop(x, pops[1], verbose = 0)),
    callrate.loc = quiet(gl.filter.callrate(x, threshold = 0.95, verbose = 0)),
    callrate.ind = quiet(gl.filter.callrate(x, method = "ind",
                                            threshold = 0.9, verbose = 0)),
    monomorphs = quiet(gl.filter.monomorphs(x, verbose = 0)),
    maf = quiet(gl.filter.maf(x, threshold = 0.05, verbose = 0)),
    rdepth = quiet(gl.filter.rdepth(x, lower = 5, upper = 50, verbose = 0)),
    reproducibility = quiet(gl.filter.reproducibility(x, threshold = 0.99,
                                                      verbose = 0)),
    secondaries = quiet(gl.filter.secondaries(x, verbose = 0)),
    allna = quiet(gl.filter.allna(x, verbose = 0)),
    subsample.loc = quiet(gl.subsample.loc(x, n = 50, replace = FALSE, verbose = 0)),
    sort = quiet(gl.sort(x, sort.by = "pop", verbose = 0))
  )
  for (nm in names(ys)) {
    expect_aligned(ys[[nm]], x)
  }
})

test_that("renaming individuals or loci cannot misalign geno.metrics", {
  x <- make_gm()
  gm <- x@other$geno.metrics
  # swapping two names: rows stay with their genotypes, by position
  y <- x
  indNames(y)[1:2] <- indNames(x)[2:1]
  locNames(y)[1:2] <- locNames(x)[2:1]
  z <- y[c(1, 3), 1:3]
  expect_identical(z@other$geno.metrics$ref.count, gm$ref.count[c(1, 3), 1:3])
  expect_null(dimnames(z@other$geno.metrics$ref.count))
  # binding does not add names either
  expect_null(dimnames(rbind(x[1:8, ], x[9:20, ])@other$geno.metrics$ref.count))
  expect_null(dimnames(cbind(x[, 1:15], x[, 16:40])@other$geno.metrics$ref.count))
})

test_that("gl.sample repeats the geno.metrics rows of individuals drawn twice", {
  x <- make_gm()
  set.seed(3)
  # more draws than individuals: rows are repeated outside "["
  capture.output(y <- gl.sample(x, nsample = 30, replace = TRUE, onepop = TRUE,
                                verbose = 0))
  expect_equal(nInd(y), 30)
  expect_aligned(y, x, ind = sub("^[0-9]+_", "", indNames(y)))
  set.seed(3)
  capture.output(y <- gl.sample(x, nsample = 10, replace = TRUE, onepop = TRUE,
                                verbose = 0))
  expect_aligned(y, x, ind = sub("^[0-9]+_", "", indNames(y)))
})
