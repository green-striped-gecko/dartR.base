# Characterization baseline for utils.check.datatype -- edge-case extension
# (function-review campaign, second-pass review of the dispatch contract).
# Complements test-utils.check.datatype.R (PR #296 baseline), which pins the
# approved diffs F1-F4 of that review; this file pins the parts of the
# contract that file does not cover: the all-all-NA courtesy-scan crash
# (PR #367's residual lead), zero-dimension and degenerate genlights,
# ploidy-vs-content classification, accept-gate matching rules, fd and
# fallback classes, FBM-backed objects, and the read-only contract.
#
# Pins CURRENT behaviour (integration-local, PR #296 merged), defects
# included. Assertions marked [approved F<n>] were flipped in Phase C of
# the v2 review (residuals branch review-utils.check.datatype-2). Where
# upstream/dev (ddaed27) behaves differently the divergence is noted
# inline:
#   - on ddaed27 the all-NA courtesy check calls gl.filter.allna, which
#     ERRORS ("Subsetting resulted in zero loci.") at verbose >= 2 when
#     every locus is all-NA, and hits "subscript out of bounds" on a
#     zero-locus object -- a function's preamble then kills the run at
#     default verbosity before the function proper sees the object.

# ---- helpers ---------------------------------------------------------------

# capture stdout in a way that survives an error inside expr
cap2 <- function(expr) {
  f <- tempfile()
  con <- file(f, "w")
  sink(con)
  err <- NA_character_
  val <- NULL
  tryCatch(val <- expr, error = function(e) err <<- conditionMessage(e))
  sink()
  close(con)
  list(lines = readLines(f, warn = FALSE), err = err, val = val)
}

# replace the genotypes of a small dartR object with matrix m (same dims)
regen <- function(g, m) {
  g@gen <- methods::new("genlight", gen = m)@gen
  adegenet::ploidy(g) <- rep(2L, nrow(m))
  g
}

# ---- the lead: all-NA courtesy scan on degenerate objects ------------------

test_that("an object with every locus all-NA survives verbose 2", {
  sub <- testset.gl[1:4, 1:6]
  m <- as.matrix(sub)
  m[] <- NA_integer_
  allna <- regen(sub, m)
  # verbose 0/1: classified from ploidy, no scan, no output
  expect_equal(utils.check.datatype(allna, verbose = 0), "SNP")
  r1 <- cap2(utils.check.datatype(allna, verbose = 1))
  expect_length(r1$lines, 0)
  expect_true(is.na(r1$err))
  # verbose 2: upstream/dev (ddaed27) ERRORS here -- its gl.filter.allna
  # courtesy call cannot return a zero-locus object ("Subsetting resulted
  # in zero loci." from the dartR '[' guard), so a caller's preamble kills
  # the run at default verbosity. PR #296's direct scan warns instead:
  r2 <- cap2(utils.check.datatype(allna, verbose = 2))
  expect_true(is.na(r2$err))
  expect_equal(r2$val, "SNP")
  expect_true(any(grepl("loci that are scored NA across all individuals",
                        r2$lines, fixed = TRUE)))
  expect_true(any(grepl("individuals that are scored NA across all loci",
                        r2$lines, fixed = TRUE)))
  # the accept gate never rescues the caller: same path with accept = "SNP"
  r3 <- cap2(utils.check.datatype(allna, accept = "SNP", verbose = 2))
  expect_true(is.na(r3$err))
})

test_that("a zero-locus (plain genlight) object survives verbose 2", {
  plain <- methods::new("genlight", gen = matrix(c(0L, 1L, 2L, 1L), nrow = 2))
  adegenet::ploidy(plain) <- c(2L, 2L)
  gl0 <- plain[, integer(0)]
  # verbose 0: classified from ploidy without touching the (empty) genotypes
  expect_equal(utils.check.datatype(gl0, verbose = 0), "SNP")
  # verbose 2: upstream/dev (ddaed27) fails with "subscript out of bounds"
  # (gl.filter.allna's 1:nLoc loop); PR #296's direct scan completes
  r <- cap2(suppressWarnings(utils.check.datatype(gl0, verbose = 2)))
  expect_true(is.na(r$err))
  expect_equal(r$val, "SNP")
})

test_that("testset.gl itself carries 3 all-NA loci that warn at verbose 2", {
  # every gl.* call on the standard fixture at default verbosity prints the
  # all-NA warning; pin the fixture fact the courtesy scan reacts to
  n_allna <- sum(colSums(!is.na(as.matrix(testset.gl))) == 0)
  expect_equal(n_allna, 3)
  r <- cap2(utils.check.datatype(testset.gl, verbose = 2))
  expect_true(any(grepl("loci that are scored NA across all individuals",
                        r$lines, fixed = TRUE)))
})

# ---- degenerate genlight objects -------------------------------------------

test_that("a genlight without ploidy fails fast with the compliance hint", {
  bare <- methods::new("genlight")
  expect_error(utils.check.datatype(bare, verbose = 0), "ploidy not set")
  # a zero-individual object has an empty ploidy slot and takes the same
  # (misleadingly worded) exit -- it is never classified
  g0 <- methods::new("genlight", gen = matrix(c(0L, 1L, 2L, 1L), nrow = 2))
  g0@gen <- list()
  g0@ploidy <- integer(0)
  g0@ind.names <- character(0)
  expect_error(utils.check.datatype(g0, verbose = 0), "ploidy not set")
})

test_that("classification is ploidy-driven, with a content consistency check", {
  # [approved F7] baseline: 0/1 (presence/absence style) values stamped
  # ploidy 2 classified as SNP -- contradictory content was not detected,
  # so a mis-stamped silico object passed an accept = 'SNP' gate into
  # dosage-based arithmetic. Now caught at the gate: uniform ploidy 2
  # with all non-missing genotypes 0/1 AND no SNP metadata (empty
  # loc.all, no SNP/SnpPosition metrics) is fatal at any verbosity. A
  # clean SNP subset lacking the homozygous-alternate class keeps its
  # loc.all slot and is unaffected (protected use case, see
  # test-gl.smearplot.R and test-gl.report.basics.R).
  m_pa <- matrix(c(0L, 1L, 1L, 0L, 0L, 1L, 1L, 0L), nrow = 2)
  g_pa <- methods::new("genlight", gen = m_pa)
  adegenet::ploidy(g_pa) <- rep(2L, 2)
  expect_error(utils.check.datatype(g_pa, accept = "SNP", verbose = 0),
               "ploidy slot and genotype content disagree")
  # [approved F7] the protected clean case: a testset.gl subset with no
  # homozygous-alternate scores passes -- its loc.all slot vouches for
  # the SNP label
  m_gl <- as.matrix(testset.gl)
  no2 <- which(colSums(m_gl == 2, na.rm = TRUE) == 0 &
                 colSums(!is.na(m_gl)) > 0)
  sub_no2 <- testset.gl[, no2[1:2]]
  expect_equal(utils.check.datatype(sub_no2, accept = "SNP", verbose = 0),
               "SNP")
  # 0/1-valued data with inferred ploidy 1 classifies as SilicoDArT
  g_bin <- methods::new("genlight", gen = matrix(c(0L, 1L, 1L, 0L), nrow = 2))
  expect_equal(utils.check.datatype(g_bin, verbose = 0), "SilicoDArT")
  # a flagless genlight built straight from a matrix (no dartR metadata at
  # all) is still classified -- flags and loc.metrics are never consulted
  g_flagless <- methods::new("genlight",
                             gen = matrix(c(0L, 1L, 2L, 1L, 0L, 2L), nrow = 2))
  expect_equal(utils.check.datatype(g_flagless, verbose = 0), "SNP")
  # [approved F7] clean SNP content (a 2 present) passes the consistency
  # check even under an accept = 'SNP' gate
  expect_equal(utils.check.datatype(g_flagless, accept = "SNP",
                                    verbose = 0), "SNP")
})

# ---- accept gate matching rules --------------------------------------------

test_that("accept matching is exact and case-sensitive", {
  expect_error(utils.check.datatype(testset.gl, accept = "snp", verbose = 0),
               "found SNP expecting snp")
  expect_error(utils.check.datatype(testset.gs, accept = "SNP", verbose = 0),
               "found SilicoDArT expecting SNP")
  expect_error(utils.check.datatype(testset.gl, accept = "SilicoDArT",
                                    verbose = 0),
               "found SNP expecting SilicoDArT")
})

# ---- non-genlight inputs ---------------------------------------------------

test_that("a data.frame is classified 'data.frame'", {
  # [approved F8] baseline: is.list() caught the data.frame before any
  # data.frame branch, so it was classified 'list' and a caller gating on
  # accept = 'matrix' read 'found list expecting matrix'. Now a data.frame
  # reports its own name.
  expect_equal(utils.check.datatype(data.frame(a = 1),
                                    accept = "data.frame",
                                    verbose = 0), "data.frame")
  expect_error(utils.check.datatype(data.frame(a = 1), accept = "list",
                                    verbose = 0), "found data.frame")
  expect_error(utils.check.datatype(data.frame(a = 1), accept = "matrix",
                                    verbose = 0), "found data.frame")
})

test_that("fd objects resolve to 'fd'; malformed fd fails fast", {
  fd <- gl.fixed.diff(testset.gl[1:10, 1:20], verbose = 0)
  expect_equal(utils.check.datatype(fd, accept = "fd", verbose = 0), "fd")
  expect_error(utils.check.datatype(fd, verbose = 0),
               "found fd")  # the default accept does not admit fd
  fd_bad <- structure(list(gl = "not a genlight"), class = "fd")
  expect_error(utils.check.datatype(fd_bad, accept = "fd", verbose = 0),
               "Fixed Difference object expected")
})

test_that("unknown classes fall through to class(x)[1] and the accept gate", {
  expect_error(utils.check.datatype(NULL, verbose = 0), "found NULL")
  expect_error(utils.check.datatype(1:3, verbose = 0), "found integer")
  # an unknown class listed in accept passes, with the warning at verbose 1
  obj <- structure(1, class = "weird")
  r0 <- cap2(utils.check.datatype(obj, accept = "weird", verbose = 0))
  expect_equal(r0$val, "weird")
  expect_length(r0$lines, 0)  # ddaed27 printed the warning even at verbose 0
  r1 <- cap2(utils.check.datatype(obj, accept = "weird", verbose = 1))
  expect_true(any(grepl("Found object of class weird", r1$lines,
                        fixed = TRUE)))
})

# ---- contracts -------------------------------------------------------------

test_that("the input object is never modified (read-only contract)", {
  before <- testset.gl
  invisible(cap2(utils.check.datatype(testset.gl, verbose = 3)))
  expect_identical(before, testset.gl)
})

test_that("verbose = 1 adds nothing over verbose = 0 (no begin/end banner)", {
  out1 <- capture.output(utils.check.datatype(testset.gs, verbose = 1))
  expect_length(out1, 0)
})

# ---- FBM-backed objects ----------------------------------------------------

test_that("an FBM-backed dartR object (empty @gen) classifies from ploidy", {
  skip_if_not(exists("gl.gen2fbm"), "gl.gen2fbm not available")
  glf <- gl.gen2fbm(testset.gl, verbose = 0)
  expect_length(glf@gen, 0)  # genotypes live in the FBM, not @gen
  expect_equal(utils.check.datatype(glf, verbose = 0), "SNP")
  # [approved F2] at verbose 2 the courtesy scan reads the FBM in column
  # blocks (no full densification) and completes
  r <- cap2(utils.check.datatype(glf, verbose = 2))
  expect_true(is.na(r$err))
  expect_equal(r$val, "SNP")
})
