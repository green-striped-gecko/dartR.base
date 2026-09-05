# Characterization tests for gl2hapmap — snapshot of CURRENT behaviour
# (review baseline, bugs included). Captured 2026-09-05 at ed99203.

test_that("gl2hapmap writes a HapMap file with the standard 11 metadata columns", {
  od <- file.path(tempdir(), "gl2hapmap_plat")
  dir.create(od, showWarnings = FALSE)
  capture.output(gl2hapmap(platypus.gl, outfile = "plat", outpath = od,
                           verbose = 0))
  hm <- read.delim(file.path(od, "plat.hmp.txt"), check.names = FALSE)

  expect_equal(nrow(hm), nLoc(platypus.gl))          # one row per locus
  expect_equal(ncol(hm), 11 + nInd(platypus.gl))     # 11 meta + 81 samples
  expect_identical(colnames(hm)[1:11],
                   c("rs#", "alleles", "chrom", "pos", "strand",
                     "assembly#", "center", "protLSID", "assayLSID",
                     "panelLSID", "QCcode"))
  expect_setequal(colnames(hm)[-(1:11)], indNames(platypus.gl))

  # alleles column keeps the A/B form; genotypes are two-letter pairs
  expect_true(all(grepl("^[ACGT]/[ACGT]$", hm$alleles)))
  gcalls <- unique(unlist(hm[, -(1:11)]))
  expect_true(all(grepl("^[ACGTN]{2}$", gcalls)))
  expect_true("NN" %in% gcalls)                      # missing code

  # platypus.gl: @chromosome is NULL -> all chrom "0".
  # Expectation updated: the baseline (captured at ed99203) emitted the
  # populated @position tag offsets verbatim; upstream dev (ddaed27, the
  # base of this branch) already zero-fills positions whose maximum is
  # < 1000 (tag offsets are not genome coordinates) -- verified on the
  # pristine base, so this diff is upstream baseline drift, not a change
  # made by this review. The zero-fill is now documented and messaged
  # (approved finding F4).
  expect_true(all(hm$chrom == 0))
  expect_true(all(hm$pos == 0))
})

test_that("gl2hapmap emits pos = 0 for every SNP when @position is NULL and pos not given", {
  od <- file.path(tempdir(), "gl2hapmap_ts")
  dir.create(od, showWarnings = FALSE)
  expect_null(testset.gl@position)
  capture.output(gl2hapmap(testset.gl, outfile = "ts", outpath = od,
                           verbose = 0))
  hm <- read.delim(file.path(od, "ts.hmp.txt"), check.names = FALSE)
  expect_true(all(hm$pos == 0))
  expect_true(all(hm$chrom == 0))
})

test_that("gl2hapmap uses nominated chrom/pos loc.metrics fields", {
  od <- file.path(tempdir(), "gl2hapmap_fields")
  dir.create(od, showWarnings = FALSE)
  capture.output(gl2hapmap(platypus.gl, outfile = "fields", outpath = od,
                           chrom = "Chrom_Platypus_Chrom_NCBIv1",
                           pos = "ChromPos_Platypus_Chrom_NCBIv1",
                           verbose = 0))
  hm <- read.delim(file.path(od, "fields.hmp.txt"), check.names = FALSE)
  # Expectation updated per approved finding F1: the nominated pos field
  # now takes precedence over the populated @position slot (previously
  # the argument was silently ignored).
  lm <- platypus.gl@other$loc.metrics
  expect_identical(sort(hm$pos),
                   sort(as.integer(lm$ChromPos_Platypus_Chrom_NCBIv1)))
  # chrom (slot NULL) comes from the nominated field
  expect_true(any(grepl("chromosome", hm$chrom)))
})

test_that("gl2hapmap pos/chrom arguments take precedence over populated slots", {
  # Expectation updated per approved finding F1: previously both
  # arguments were silently ignored whenever the corresponding slot was
  # populated (pos stayed 77, chrom stayed 'Z').
  od <- file.path(tempdir(), "gl2hapmap_ignored")
  dir.create(od, showWarnings = FALSE)
  p2 <- platypus.gl
  p2@position <- rep(77L, nLoc(p2))
  p2@chromosome <- factor(rep("Z", nLoc(p2)))
  capture.output(gl2hapmap(p2, outfile = "ign", outpath = od,
                           chrom = "Chrom_Platypus_Chrom_NCBIv1",
                           pos = "SnpPosition", verbose = 0))
  hm <- read.delim(file.path(od, "ign.hmp.txt"), check.names = FALSE)
  lm <- platypus.gl@other$loc.metrics
  expect_identical(sort(hm$pos), sort(as.integer(lm$SnpPosition)))
  expect_false(any(hm$chrom == "Z"))
  expect_true(any(grepl("chromosome", hm$chrom)))
})

test_that("gl2hapmap rejects SilicoDArT data with an informative condition", {
  # Expectation updated per approved finding F6: the former
  # cat(error()); stop() raised a condition with an empty message; the
  # rejection now goes through utils.check.datatype(accept = 'SNP').
  od <- file.path(tempdir(), "gl2hapmap_gs")
  dir.create(od, showWarnings = FALSE)
  err <- tryCatch(capture.output(gl2hapmap(testset.gs, outpath = od,
                                           verbose = 0)),
                  error = function(e) e)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "SilicoDArT")
})

test_that("gl2hapmap is silent at verbose = 0", {
  # Expectation updated per approved findings F2 (saved-file message now
  # gated at verbose >= 2) and F5 (invisible NULL return).
  od <- file.path(tempdir(), "gl2hapmap_v0")
  dir.create(od, showWarnings = FALSE)
  out <- capture.output(gl2hapmap(platypus.gl, outfile = "v0", outpath = od,
                                  verbose = 0))
  expect_length(out, 0)
})

test_that("gl2hapmap fails fast with a clear message on a genlight without loc.all", {
  # Expectation updated per approved finding F3: previously died deep in
  # the recode with 'number of items to replace is not a multiple of
  # replacement length'.
  gl.sim <- adegenet::glSim(10, 20, ploidy = 2)
  od <- file.path(tempdir(), "gl2hapmap_sim")
  dir.create(od, showWarnings = FALSE)
  expect_error(capture.output(gl2hapmap(gl.sim, outfile = "sim",
                                        outpath = od, verbose = 0)),
               "allele definitions")
})
