# Characterization tests for gl2gi — snapshot of CURRENT behaviour
# (review baseline, bugs included). Captured 2026-09-05 at ed99203.

test_that("gl2gi converts SNP genlight to genind preserving names, pop, ploidy", {
  gi <- suppressWarnings(gl2gi(testset.gl, verbose = 0))

  expect_s4_class(gi, "genind")
  expect_equal(unique(adegenet::ploidy(gi)), 2)
  expect_equal(adegenet::nInd(gi), 274)
  expect_identical(adegenet::indNames(gi), indNames(testset.gl))
  expect_identical(as.character(adegenet::pop(gi)),
                   as.character(pop(testset.gl)))

  # Expectation updated per approved finding F1: all-NA loci are now
  # removed up front and loc.metrics is subset to the surviving loci, so
  # the 752-locus genind carries 752-row locus metadata (was 755).
  expect_equal(adegenet::nLoc(gi), 752)
  expect_equal(nrow(gi@other$loc.metrics), 752)

  # surviving locus names are unmangled and in order
  kept <- locNames(testset.gl)[locNames(testset.gl) %in%
                                 adegenet::locNames(gi)]
  expect_identical(adegenet::locNames(gi), kept)

  # other slot survives wholesale
  expect_true(all(c("loc.metrics", "ind.metrics", "latlon") %in%
                    names(gi@other)))
  expect_equal(nrow(gi@other$ind.metrics), 274)
})

test_that("gl2gi genotype coding: 0 -> hom ref allele, 1 -> het, 2 -> hom alt allele", {
  x <- gl.filter.allna(testset.gl, verbose = 0)[1:20, 1:50]
  gi <- suppressWarnings(gl2gi(x, verbose = 0))
  m <- as.matrix(x)
  tab <- adegenet::tab(gi)

  # loci monomorphic within the subset carry only one allele column in
  # the genind; check the coding on loci where both alleles survive.
  # For those, the first allele of loc.all is counted 2/1/0 for
  # genotypes 0/1/2 (dartR 0 = homozygous reference), and NA stays NA.
  checked <- 0L
  for (j in seq_len(nLoc(x))) {
    a1 <- substr(x@loc.all[j], 1, 1)
    a2 <- substr(x@loc.all[j], 3, 3)
    c1 <- paste0(locNames(x)[j], ".", a1)
    c2 <- paste0(locNames(x)[j], ".", a2)
    if (!(c1 %in% colnames(tab)) || !(c2 %in% colnames(tab))) next
    k <- !is.na(m[, j])
    expect_equal(unname(tab[k, c1]), unname(2 - m[k, j]))
    expect_true(all(is.na(tab[!k, c1])))
    checked <- checked + 1L
    if (checked >= 5L) break
  }
  expect_gte(checked, 3L)
})

test_that("gl2gi rejects SilicoDArT input", {
  # Expectation updated per approved finding F2 (DAT7): accept = 'SNP'
  # now rejects presence/absence data, which previously converted to a
  # meaningless ploidy-2 genind.
  expect_error(gl2gi(testset.gs, verbose = 0))
})

test_that("gl2gi fabricates allele labels when loc.all is NULL (current behaviour)", {
  gl.sim <- adegenet::glSim(8, 15, ploidy = 2)
  gi <- suppressWarnings(gl2gi(gl.sim, verbose = 0))
  alleles <- unlist(lapply(gi@all.names, paste, collapse = "/"))
  # locus 1 is labelled C/G, the rest A/T (or a subset when only one
  # allele is observed) — invented, not read from the data
  expect_true(all(unlist(gi@all.names) %in% c("A", "T", "C", "G")))
})

test_that("gl2gi converts a single-locus genlight", {
  # Expectation updated per approved finding F3: the drop = FALSE guard
  # fixes the former crash ("X is not a matrix") on one-locus objects.
  # Note: df2genind removes the one individual with no scored genotype
  # at this locus (274 -> 273); its ind.metrics row is not removed —
  # pre-existing behaviour, recorded in the report as a follow-up.
  x <- gl.filter.allna(testset.gl, verbose = 0)[, 1]
  gi <- suppressWarnings(gl2gi(x, verbose = 0))
  expect_s4_class(gi, "genind")
  expect_equal(adegenet::nLoc(gi), 1)
  expect_equal(adegenet::nInd(gi), 273)
})

test_that("gl2gi is silent at verbose = 0 (console output)", {
  out <- capture.output(gi <- suppressWarnings(gl2gi(testset.gl,
                                                     verbose = 0)))
  expect_length(out, 0)
})
