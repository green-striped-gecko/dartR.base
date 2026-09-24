# Characterization tests for rbind.dartR (R/utils.dartR.class.def.r)
# Baseline snapshotted before review (origin/dev at 1dd6a3a). Assertions
# marked [approved diff N] were flipped in Phase C against change N in
# function-review/reports/dartR.base/rbind.dartR.md.

test_that("same loci, same order: genotypes and dimensions", {
  x <- platypus.gl
  r <- rbind(x[1:5, ], x[6:10, ])
  expect_s4_class(r, "dartR")
  expect_equal(nInd(r), 10L)
  expect_equal(nLoc(r), nLoc(x))
  expect_identical(as.matrix(r), as.matrix(x[1:10, ]))
  expect_identical(indNames(r), indNames(x)[1:10])
})

test_that("metadata is carried over [approved diff 2]", {
  x <- platypus.gl
  r <- rbind(x[1:5, ], x[6:10, ])
  # baseline: @other was empty on the in-memory path
  expect_equal(nrow(r@other$loc.metrics), nLoc(x))
  expect_identical(r@other$loc.metrics$AlleleID, x@other$loc.metrics$AlleleID)
  expect_equal(nrow(r@other$ind.metrics), 10L)
  expect_identical(as.character(r@other$ind.metrics$id), indNames(r))
  expect_equal(nrow(r@other$latlon), 10L)
  expect_false(any(unlist(r@other$loc.metrics.flags)))
  h <- r@other$history
  expect_equal(length(h), length(x[1:5, ]@other$history) + 1L)
  expect_identical(deparse(h[[length(h)]]), "rbind.dartR(n.objects = 2)")
})

test_that("ind.metrics columns missing from one object are NA", {
  x <- platypus.gl
  a <- x[1:3, ]
  b <- x[4:6, ]
  a@other$ind.metrics$extra <- 1:3
  r <- rbind(a, b)
  expect_equal(r@other$ind.metrics$extra, c(1:3, NA, NA, NA))
})

test_that("loci in another order are aligned by name [approved diff 1]", {
  x <- platypus.gl
  set.seed(1)
  ord <- sample(nLoc(x))
  r <- rbind(x[1:5, ], x[6:10, ord])
  expect_identical(locNames(r), locNames(x))
  # baseline: 1722 of 4636 genotypes of rows 6-10 were misplaced
  expect_identical(as.matrix(r), as.matrix(x[1:10, ]))
})

test_that("different allele coding stops [approved diff 3]", {
  x <- platypus.gl
  b <- x[6:10, ]
  b@loc.all <- vapply(strsplit(b@loc.all, "/"),
                      function(a) paste(rev(a), collapse = "/"), "")
  expect_error(rbind(x[1:5, ], b), "codes the alleles")
})

test_that("SNP with SilicoDArT stops [approved diff 4]", {
  x <- platypus.gl[1:3, 1:255]
  gs <- testset.gs[1:3, ]
  locNames(gs) <- locNames(x)
  # baseline: combined silently into one object with ploidy 2 and 1
  expect_error(rbind(x, gs), "mix ploidy")
})

test_that("arguments: NULL allowed, others stop [approved diff 5]", {
  x <- platypus.gl
  expect_equal(nInd(rbind(NULL, x[1:2, ])), 2L)
  expect_equal(nInd(do.call(rbind, list(x[1:2, ], x[3:4, ]))), 4L)
  # baseline: a non-genlight argument was silently ignored
  expect_error(rbind(x[1:2, ], data.frame(a = 1)), "accepts only")
})

test_that("FBM path gives the same genotypes and metadata", {
  skip_if_not_installed("bigstatsr")
  x <- platypus.gl
  r <- rbind(gl.gen2fbm(x[1:5, ], verbose = 0),
             gl.gen2fbm(x[6:10, ], verbose = 0))
  expect_equal(unname(as.matrix(r)), unname(as.matrix(x[1:10, ])))
  # baseline: ind.metrics kept the first object's 5 rows for 10 individuals
  expect_equal(nrow(r@other$ind.metrics), 10L)
  expect_false(any(unlist(r@other$loc.metrics.flags)))
})
