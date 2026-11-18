# Characterization baseline (function-review, 2026-09-24): outputs of the
# reviewed state. Detects change; does not assert correctness.
test_that("baseline: gen -> FBM -> gen round trip preserves the object", {
  skip_if_not_installed("bigstatsr")
  x <- testset.gl
  capture.output(f <- gl.gen2fbm(x, verbose = 0))
  expect_length(f@gen, 0)
  expect_s4_class(f@fbm, "FBM.code256")
  expect_equal(unname(as.matrix(f)) + 0, unname(as.matrix(x)) + 0)
  capture.output(b <- gl.fbm2gen(f, verbose = 0))
  expect_null(b@fbm)
  expect_identical(as.matrix(b), as.matrix(x))
  for (s in c("ind.names", "loc.names", "loc.all", "ploidy", "pop",
              "position", "chromosome", "n.loc", "other")) {
    expect_identical(slot(b, s), slot(x, s), label = s)
  }
  # small blocks give the same result
  capture.output(f7 <- gl.gen2fbm(x, chunk = 7L, verbose = 0))
  expect_equal(unname(as.matrix(f7)) + 0, unname(as.matrix(x)) + 0)
  # already FBM-backed: returned as is
  capture.output(f2 <- gl.gen2fbm(f, verbose = 0))
  expect_identical(f2@fbm$backingfile, f@fbm$backingfile)
})

test_that("baseline: unsupported inputs", {
  skip_if_not_installed("bigstatsr")
  expect_error(capture.output(gl.gen2fbm(testset.gs, verbose = 0)),
               "Only SNP data")
  # [approved diff, change 2] a plain genlight failed in both functions
  g <- new("genlight", as.matrix(testset.gl)[1:20, 1:30], ploidy = 2)
  capture.output(gf <- gl.gen2fbm(g, verbose = 0))
  expect_s4_class(gf, "dartR")
  expect_equal(unname(as.matrix(gf)) + 0, unname(as.matrix(g)) + 0)
  capture.output(gb <- gl.fbm2gen(g, verbose = 0))
  expect_identical(gb, g)
  # gen-backed input to fbm2gen is returned unchanged
  capture.output(b <- gl.fbm2gen(testset.gl, verbose = 0))
  expect_identical(b, testset.gl)
})

test_that("gl.fbm2gen gives the same object for any block size", {
  skip_if_not_installed("bigstatsr")
  x <- platypus.gl
  capture.output(f <- gl.gen2fbm(x, verbose = 0))
  ref <- new("genlight", as.matrix(x), ploidy = x@ploidy)@gen
  for (k in c(8L, 16L, 504L, 1000L, 1008L)) {
    gen <- dartR.base:::utils.fbm2gen.copy(f@fbm, x@ploidy, loc_per = k)
    expect_identical(gen, ref, label = paste("loci per block", k))
  }
  capture.output(b <- gl.fbm2gen(f, verbose = 0))
  expect_identical(as.matrix(b), as.matrix(x))
  expect_identical(b@ploidy, x@ploidy)
  expect_null(b@fbm)
  # deprecated 'chunk' is accepted, ignored and reported
  capture.output(b2 <- gl.fbm2gen(f, chunk = 7L, verbose = 0))
  expect_identical(b2@gen, b@gen)
  expect_output(gl.fbm2gen(f, chunk = 7L, verbose = 1), "deprecated")
})

test_that("gl.fbm2gen builds the same SNPbins as adegenet", {
  # The SNPbins are packed directly from the FBM; they must be identical to
  # those of adegenet's constructor, including its rules for the number of
  # bit planes and the SNPbin ploidy (an individual without a genotype 2
  # gets one plane and ploidy 1; an all-missing individual ploidy NA)
  same <- function(x) {
    capture.output(y <- gl.fbm2gen(gl.gen2fbm(x, verbose = 0), verbose = 0))
    ref <- new("genlight", as.matrix(x), ploidy = x@ploidy)@gen
    expect_identical(y@gen, ref)
    expect_identical(as.matrix(y), as.matrix(x))
    expect_identical(y@ploidy, x@ploidy)
  }
  mk <- function(n, p, ploidy = 2L, na = 0.1) {
    m <- matrix(sample(0:2, n * p, TRUE), n)
    m[ploidy == 1L, ] <- m[ploidy == 1L, ] %% 2L
    m[sample(length(m), floor(length(m) * na))] <- NA
    new("genlight", m, ploidy = ploidy)
  }
  set.seed(5)
  same(testset.gl)
  for (p in c(1, 7, 8, 9, 17, 100)) {
    same(mk(1, p))
    same(mk(3, p))
  }
  same(mk(5, 50, ploidy = c(1L, 2L, 1L, 2L, 2L)))
  same(mk(5, 50, na = 0))
  same(mk(5, 50, na = 0.9))
  m <- rbind(c(NA, NA, NA), c(0L, 0L, NA), c(1L, 0L, 1L), c(2L, 1L, 0L))
  same(new("genlight", m, ploidy = 2L))
  same(new("genlight", m[3:4, ], ploidy = c(1L, 2L)))
  # stored SNPbins with extra padding bytes (as in platypus.gl) decode the
  # same; the result has adegenet's standard layout
  same(platypus.gl)
})

test_that("gl.fbm2gen refuses values that SNPbins cannot hold", {
  g <- new("genlight", rbind(c(2L, 0L, 1L), c(0L, 1L, NA)), ploidy = 2L)
  capture.output(f <- gl.gen2fbm(g, verbose = 0))
  f@ploidy[1] <- 1L                      # a haploid individual with a 2
  expect_error(gl.fbm2gen(f, verbose = 0), "haploid")
  G <- bigstatsr::FBM(2, 3, init = 0.5)
  expect_error(dartR.base:::utils.fbm2gen.copy(G, 2L, loc_per = 8L),
               "other than 0, 1, 2")
})

test_that("an object with both slots populated is reported as invalid", {
  skip_if_not_installed("bigstatsr")
  capture.output(f <- gl.gen2fbm(testset.gl, verbose = 0))
  bad <- f
  bad@gen <- testset.gl@gen
  expect_error(capture.output(gl.gen2fbm(bad, verbose = 0)),
               "both @fbm and @gen")
})

test_that("polyploid dosages are refused, not silently set to missing", {
  tg <- new("genlight",
            matrix(c(0, 1, 2, 3, 4, 1, 2, 0, 1, 4, 3, NA), nrow = 2,
                   byrow = TRUE),
            ploidy = 4L)
  before <- list.files(tempdir(), pattern = "^geno_")
  expect_error(gl.gen2fbm(tg, verbose = 0), "diploid|above 2")
  # no new file left behind (older ones may be cleaned up meanwhile)
  expect_length(setdiff(list.files(tempdir(), pattern = "^geno_"), before), 0)
})

test_that("gen -> FBM copies the decoded genotypes exactly (SNPbin layout)", {
  # The copy reads the packed SNPbin bits directly; compare with adegenet's
  # own decoding so that a change of the SNPbin layout is caught.
  same <- function(x, ...) {
    capture.output(f <- gl.gen2fbm(x, verbose = 0, ...))
    m <- as.matrix(x)
    expect_identical(f@fbm[seq_len(nrow(m)), seq_len(ncol(m)), drop = FALSE],
                     unname(m) + 0)
  }
  mk <- function(n, p, ploidy = 2L, na = 0.1) {
    m <- matrix(sample(0:2, n * p, TRUE), n)
    m[ploidy == 1L, ] <- m[ploidy == 1L, ] %% 2L
    m[sample(length(m), floor(length(m) * na))] <- NA
    new("genlight", m, ploidy = ploidy)
  }
  set.seed(11)
  same(testset.gl)
  # loci counts around the 8 loci per byte, one individual, one locus
  for (p in c(1, 7, 8, 9, 17, 100)) {
    same(mk(1, p))
    same(mk(3, p))
  }
  same(mk(5, 50, ploidy = c(1L, 2L, 1L, 2L, 2L)))   # haploid individuals
  same(mk(5, 50, na = 0))                           # no missing data
  same(mk(5, 50, na = 0.9))                         # mostly missing
  m0 <- matrix(NA_integer_, 3, 10); m0[, 1] <- 0L    # all-missing loci
  same(new("genlight", m0, ploidy = 2L))
  # block boundaries: tiny blocks and the deprecated 'chunk'
  for (k in c(1L, 8L, 9L, 64L)) same(testset.gl, chunk = k)
  op <- options(bigstatsr.block.sizeGB = 1e-5)
  on.exit(options(op))
  expect_lt(bigstatsr::block_size(nInd(testset.gl)), nLoc(testset.gl))
  same(testset.gl)
})

test_that("the deprecated 'chunk' argument is reported", {
  expect_output(gl.gen2fbm(testset.gl, chunk = 100L, verbose = 1),
                "deprecated")
})
