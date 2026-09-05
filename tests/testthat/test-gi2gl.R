# Characterization tests for gi2gl — snapshot of CURRENT behaviour
# (review baseline, bugs included), including the gl2gi -> gi2gl round
# trip. Captured 2026-09-05 at ed99203.

test_that("round trip gl -> gi -> gl on raw testset2.gl", {
  # Baseline behaviour: crashed via gl2gi's loc.metrics desync on dropped
  # all-NA loci (gl2gi finding F1, fixed on the separate review-gl2gi
  # branch). Written to hold in both merge states: with the desync still
  # present gi2gl errors on the row-count mismatch; with gl2gi fixed the
  # genind is clean and converts.
  gi <- suppressWarnings(gl2gi(testset2.gl, verbose = 0))
  if (nrow(gi@other$loc.metrics) == adegenet::nLoc(gi)) {
    rt <- gi2gl(gi, verbose = 0)
    expect_equal(nLoc(rt), adegenet::nLoc(gi))
  } else {
    expect_error(gi2gl(gi, verbose = 0), "752|755")
  }
})

test_that("round trip on clean data preserves structure; dosage polarity flips per locus", {
  x <- gl.filter.allna(testset2.gl, verbose = 0)   # 752 loci, no all-NA
  gi <- suppressWarnings(gl2gi(x, verbose = 0))
  rt <- gi2gl(gi, verbose = 0)

  # structure and identity survive
  expect_true(inherits(rt, "genlight"))
  expect_equal(unique(ploidy(rt)), 2)
  expect_equal(nInd(rt), nInd(x))
  expect_equal(nLoc(rt), nLoc(x))
  expect_identical(indNames(rt), indNames(x))
  expect_identical(locNames(rt), locNames(x))
  expect_identical(as.character(pop(rt)), as.character(pop(x)))
  expect_equal(nrow(rt@other$loc.metrics), 752)
  expect_equal(nrow(rt@other$ind.metrics), 274)

  # NA pattern survives exactly
  m0 <- unname(as.matrix(x))
  m2 <- unname(as.matrix(rt))
  expect_identical(is.na(m0), is.na(m2))

  # genotypes do NOT survive verbatim: gi2gl counts the first allele in
  # the genind tab, so loci whose original alt allele sorts first come
  # back complemented (0 <-> 2, hets unchanged). On this data 596 of 752
  # loci flip; every locus is either identical or an exact 2-x
  # complement (the documented reference-allele ambiguity).
  flipped <- 0L
  for (j in seq_len(ncol(m0))) {
    k <- !is.na(m0[, j])
    if (all(m2[k, j] == m0[k, j])) next
    expect_equal(m2[k, j], 2 - m0[k, j])
    flipped <- flipped + 1L
  }
  expect_equal(flipped, 596L)

  # Expectation updated per approved finding F4: loc.all is now
  # reconstructed from the genind allele names instead of a uniform
  # placeholder. For loci biallelic in the data the pair matches the
  # original up to order (flipped loci carry the reversed spelling,
  # consistent with their complemented dosage); loci monomorphic in the
  # data come back as 'a/a'. Not identical overall because of the
  # reversals.
  expect_false(identical(rt@loc.all, x@loc.all))
  n2 <- lengths(gi@all.names) == 2
  same_set <- mapply(function(o, r) setequal(o, r),
                     strsplit(x@loc.all, "/"),
                     strsplit(rt@loc.all, "/"))
  expect_true(all(same_set[n2]))

  # history is appended on the returned object
  expect_gt(length(rt@other$history), length(x@other$history))
})

test_that("gi2gl rejects non-genind input", {
  expect_error(gi2gl(testset2.gl, verbose = 0), "genind object required")
})

test_that("gi2gl converts a single-locus genind", {
  # Expectation updated per approved finding F2: the cumsum offset walk
  # fixes the former crash ("argument is of length zero") on one-locus
  # genind objects.
  df1 <- data.frame(L1 = c("A/A", "A/T", "T/T", "A/A"),
                    row.names = paste0("i", 1:4))
  gi1 <- adegenet::df2genind(df1, sep = "/")
  g1 <- gi2gl(gi1, verbose = 0)
  expect_true(inherits(g1, "genlight"))
  expect_equal(nLoc(g1), 1)
  expect_equal(nInd(g1), 4)
})

test_that("gi2gl rejects a genind with more than two alleles at a locus", {
  # Expectation updated per approved finding F1: previously the
  # column-offset walk assumed <= 2 alleles and every locus after the
  # multiallelic one silently read the wrong tab column (loc2 came back
  # as loc1's third-allele counts, 0,0,1,2). Now a fatal error names the
  # offending locus.
  df <- data.frame(loc1 = c("101/101", "101/102", "102/103", "103/103"),
                   loc2 = c("201/201", "201/202", "202/202", "201/201"),
                   row.names = paste0("ind", 1:4))
  gi.ms <- adegenet::df2genind(df, sep = "/")
  expect_error(gi2gl(gi.ms, verbose = 0), "loc1")
})

test_that("gi2gl is silent at verbose = 0", {
  x <- gl.filter.allna(testset2.gl, verbose = 0)
  gi <- suppressWarnings(gl2gi(x, verbose = 0))
  out <- capture.output(rt <- gi2gl(gi, verbose = 0))
  expect_length(out, 0)
})
