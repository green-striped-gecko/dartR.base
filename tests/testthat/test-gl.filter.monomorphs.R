test_that("gl.filter.monomorphs removes exactly the monomorphic + all-NA loci (SNP)", {
  x <- testset.gl
  m <- as.matrix(x)
  is_mono <- apply(m, 2, function(col) {
    v <- col[!is.na(col)]
    length(v) == 0 || all(v == 0) || all(v == 2)
  })
  x2 <- gl.filter.monomorphs(x, verbose = 0)
  expect_equal(nLoc(x2), nLoc(x) - sum(is_mono))
  expect_setequal(locNames(x2), locNames(x)[!is_mono])
  expect_equal(nrow(x2@other$loc.metrics), nLoc(x2))
  expect_equal(nInd(x2), nInd(x))
  expect_true(all(ploidy(x2) == 2))
  expect_true(x2@other$loc.metrics.flags$monomorphs)
  expect_length(x2@other$history, length(x@other$history) + 1)
})

test_that("gl.filter.monomorphs removes exactly the monomorphic + all-NA loci (SilicoDArT)", {
  x <- testset.gs
  m <- as.matrix(x)
  is_mono <- apply(m, 2, function(col) {
    v <- col[!is.na(col)]
    length(v) == 0 || all(v == 0) || all(v == 1)
  })
  x2 <- gl.filter.monomorphs(x, verbose = 0)
  expect_equal(nLoc(x2), nLoc(x) - sum(is_mono))
  expect_setequal(locNames(x2), locNames(x)[!is_mono])
  expect_true(all(ploidy(x2) == 1))
})

test_that("all-heterozygous loci (uniform dosage 1) are correctly retained as polymorphic", {
  m2 <- matrix(c(1, 0, 1, 1, 2, 1, 1, 1, 1, 1, 0, 1), nrow = 3, byrow = TRUE,
               dimnames = list(paste0("I", 1:3), paste0("L", 1:4)))
  # L1: 1,1,1 all het -> polymorphic (retain); L2: 0,2,1 -> polymorphic;
  # L3: 1,1,0 -> polymorphic; L4: 1,1,1 all het -> retain
  gl_small <- new("genlight", gen = m2, ploidy = 2)
  pop(gl_small) <- factor(rep("p1", 3))
  gl_small <- gl.compliance.check(gl_small, verbose = 0)
  x2 <- gl.filter.monomorphs(gl_small, verbose = 0)
  expect_equal(nLoc(x2), 4)
})

test_that("a run with no monomorphs present returns the object with flag TRUE", {
  x <- gl.filter.monomorphs(testset.gl, verbose = 0)
  x2 <- gl.filter.monomorphs(x, verbose = 0)
  expect_equal(nLoc(x2), nLoc(x))
  expect_true(x2@other$loc.metrics.flags$monomorphs)
})

test_that("verbose = 0 is fully silent, including on an unassigned call (F1 fix)", {
  x <- testset.gl
  out0 <- capture.output(gl.filter.monomorphs(x, verbose = 0))
  expect_length(out0, 0)
})

test_that("history does not carry the internal gl.drop.loc entry (F2 fix)", {
  x <- testset.gl
  x2 <- gl.filter.monomorphs(x, verbose = 0)
  new_entries <- x2@other$history[(length(x@other$history) + 1):length(x2@other$history)]
  expect_length(new_entries, 1)
  expect_true(grepl("gl.filter.monomorphs", deparse(new_entries[[1]])[1]))
  expect_false(any(grepl("gl.drop.loc", vapply(new_entries, function(h) deparse(h)[1], ""))))
})

# ---------------------------------------------------------------------------
# File-backed (FBM) objects (2026-11-10). Monomorphic loci are found from
# genotype counts read from the backing file in blocks
# (bigstatsr::big_counts), where as.matrix(x) read the whole matrix into
# memory. The object returned must be the one the in-memory code returns.

fbm_mono_case <- function() {
  set.seed(3)
  m <- matrix(sample(0:2, 12 * 30, replace = TRUE), 12, 30,
              dimnames = list(sprintf("ind%02d", 1:12), sprintf("L%02d", 1:30)))
  m[sample.int(length(m), 40)] <- NA
  m[, 3] <- NA             # all missing
  m[, 9] <- 0              # fixed for the reference allele
  m[, 10] <- 2             # fixed for the alternate allele
  m[, 11] <- 1             # all heterozygous: polymorphic
  m[, 12] <- c(0, NA)      # 0 and missing: monomorphic
  m[, 13] <- c(2, 2, NA)   # 2 and missing: monomorphic
  m[, 14] <- c(0, 2)       # both homozygotes: polymorphic
  x <- new("genlight", gen = m, ploidy = 2)
  pop(x) <- factor(rep(c("A", "B", "C"), each = 4))
  gl.compliance.check(x, verbose = 0)
}

test_that("FBM-backed object: same loci removed as in memory (constructed cases)", {
  skip_if_not_installed("bigstatsr")
  x <- fbm_mono_case()
  fb <- gl.gen2fbm(x, verbose = 0)
  a <- gl.filter.monomorphs(x, verbose = 0)
  b <- gl.filter.monomorphs(fb, verbose = 0)
  expect_setequal(setdiff(locNames(x), locNames(b)),
                  c("L03", "L09", "L10", "L12", "L13"))
  expect_identical(locNames(b), locNames(a))
  expect_identical(indNames(b), indNames(a))
  expect_equal(unname(as.matrix(b)), unname(as.matrix(a)))
  expect_identical(b@other$loc.metrics, a@other$loc.metrics)
  # still file-backed, flag set, one history entry
  expect_false(is.null(b@fbm))
  expect_true(b@other$loc.metrics.flags$monomorphs)
  expect_length(b@other$history, length(fb@other$history) + 1)
})

test_that("FBM-backed object: same loci removed and same messages as in memory (testset.gl)", {
  skip_if_not_installed("bigstatsr")
  fb <- gl.gen2fbm(testset.gl, verbose = 0)
  out.a <- capture.output(a <- gl.filter.monomorphs(testset.gl, verbose = 3))
  out.b <- capture.output(b <- gl.filter.monomorphs(fb, verbose = 3))
  expect_identical(locNames(b), locNames(a))
  expect_equal(unname(as.matrix(b)), unname(as.matrix(a)))
  # counts of monomorphic and all-NA loci are reported the same
  expect_identical(grep("loci|individuals|populations", out.b, value = TRUE),
                   grep("loci|individuals|populations", out.a, value = TRUE))
})
