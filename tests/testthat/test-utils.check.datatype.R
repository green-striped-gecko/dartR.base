# Characterization tests for utils.check.datatype
# Baseline snapshotted before review (review-utils.check.datatype).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("dispatch across classes", {
  expect_equal(utils.check.datatype(testset2.gl, verbose = 0), "SNP")
  expect_equal(utils.check.datatype(testset2.gs, verbose = 0), "SilicoDArT")
  expect_equal(utils.check.datatype(dist(matrix(1:9, 3)),
        accept = "dist", verbose = 0), "dist")
  expect_equal(utils.check.datatype(matrix(1:4, 2),
        accept = "matrix", verbose = 0), "matrix")
  expect_equal(utils.check.datatype(list(a = 1),
        accept = "list", verbose = 0), "list")
  v <- withVisible(utils.check.datatype(testset2.gl, verbose = 0))
  expect_false(v$visible)
  expect_error(utils.check.datatype(testset2.gl, accept = "dist",
                                    verbose = 0))
})

test_that("abnormal ploidy handling", {
  gg <- testset2.gl[1:10, 1:20]
  gg@ploidy <- as.integer(c(rep(2, 5), rep(1, 5)))
  # classified SNP (polyploid paths depend on non-fatal handling)
  expect_equal(suppressWarnings(utils.check.datatype(gg, verbose = 0)),
               "SNP")
  o <- capture.output(r <- utils.check.datatype(gg, verbose = 2))
  # [approved diff F2] baseline: no notice of the abnormal ploidy at
  # any verbosity. Now announced at verbose >= 2.
  expect_true(any(grepl("ploidy is not uniformly 2", o)))
  # and still silent at verbose 0
  o0 <- capture.output(r0 <- utils.check.datatype(gg, verbose = 0))
  expect_length(o0, 0)
})

test_that("all-NA warnings at verbose 2", {
  # all-NA individual
  set.seed(7)
  mm <- matrix(sample(0:2, 200, replace = TRUE), nrow = 10)
  mm[3, ] <- NA
  g3 <- new("genlight", gen = mm, ploidy = 2)
  indNames(g3) <- paste0("i", 1:10); locNames(g3) <- paste0("l", 1:20)
  pop(g3) <- factor(rep("A", 10))
  g3 <- suppressWarnings(gl.compliance.check(g3, verbose = 0))
  o <- capture.output(utils.check.datatype(g3, verbose = 2))
  # [approved diff F1] baseline: the individuals warning carried the
  # LOCI wording.
  expect_false(any(grepl("loci that are scored NA", o)))
  expect_true(any(grepl("individuals that are scored NA across all loci", o)))
  # all-NA loci (testset2.gl has 3): warning present at verbose 2
  o2 <- capture.output(utils.check.datatype(testset2.gl, verbose = 2))
  expect_true(any(grepl("loci that are scored NA", o2)))
  # and silent at verbose 0
  o0 <- capture.output(utils.check.datatype(testset2.gl, verbose = 0))
  expect_length(o0, 0)
})

test_that("unknown class warning", {
  # [approved diff F4] baseline: warned at verbose 0 before the
  # accept-check fatal. Now gated at verbose >= 1.
  o <- capture.output(r <- try(utils.check.datatype("a string",
        verbose = 0), silent = TRUE))
  expect_length(o, 0)
  expect_true(inherits(r, "try-error"))
})

test_that("accept = 'genlight' semantics", {
  # [approved diff F3] baseline: fatal - datatype is never literally
  # 'genlight', so a genlight-only accept rejected every genlight
  # object. Now 'genlight'/'dartR' admit both genotype datatypes,
  # unless a specific datatype is also listed (which then governs).
  expect_equal(utils.check.datatype(testset2.gl, accept = "genlight",
                                    verbose = 0), "SNP")
  expect_equal(utils.check.datatype(testset2.gs, accept = "dartR",
                                    verbose = 0), "SilicoDArT")
  # the specific listing still governs: c('genlight','SNP') is SNP-only
  expect_error(utils.check.datatype(testset2.gs,
        accept = c("genlight", "SNP"), verbose = 0))
})
