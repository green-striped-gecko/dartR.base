# Characterization tests for gl2geno — snapshot of CURRENT behaviour
# (review baseline, bugs included). Captured 2026-09-05 at ed99203.

test_that("gl2geno writes geno and lfmm files for SNP data (all-NA loci dropped)", {
  od <- file.path(tempdir(), "gl2geno_snp")
  dir.create(od, showWarnings = FALSE)
  res <- gl2geno(testset2.gl, outfile = "ts", outpath = od, verbose = 0)

  # returns NULL invisibly
  expect_null(res)

  geno <- readLines(file.path(od, "ts.geno"))
  lfmm <- readLines(file.path(od, "ts.lfmm"))

  # testset2.gl has 755 loci of which 3 are all-NA; gl2geno removes them
  # via gl.filter.allna before writing
  expect_equal(length(geno), 752)          # geno: one row per locus
  expect_equal(length(lfmm), 274)          # lfmm: one row per individual
  expect_equal(nchar(geno[1]), 274)        # geno rows: one char per individual

  # geno alphabet is 0/1/2 with 9 as the missing code, no separators
  chars <- sort(unique(strsplit(paste(geno, collapse = ""), "")[[1]]))
  expect_equal(chars, c("0", "1", "2", "9"))

  # lfmm is space-separated, same alphabet
  f1 <- strsplit(lfmm[1], " ")[[1]]
  expect_equal(length(f1), 752)
  expect_true(all(f1 %in% c("0", "1", "2", "9")))

  # values agree with the genotype matrix (NA -> 9), individuals in order
  m <- as.matrix(gl.filter.allna(testset2.gl, verbose = 0))
  m[is.na(m)] <- 9
  expect_equal(as.integer(f1), unname(m[1, ]))
})

test_that("gl2geno handles SilicoDArT data (0/1 alphabet plus 9)", {
  od <- file.path(tempdir(), "gl2geno_gs")
  dir.create(od, showWarnings = FALSE)
  gl2geno(testset2.gs, outfile = "gs", outpath = od, verbose = 0)
  geno <- readLines(file.path(od, "gs.geno"))
  expect_equal(length(geno), 752)          # testset2.gs also has 3 all-NA loci
  chars <- sort(unique(strsplit(paste(geno, collapse = ""), "")[[1]]))
  expect_equal(chars, c("0", "1", "9"))
})

test_that("gl2geno is silent at verbose = 0", {
  od <- file.path(tempdir(), "gl2geno_v0")
  dir.create(od, showWarnings = FALSE)
  out <- capture.output(gl2geno(testset2.gl, outfile = "v0", outpath = od,
                                verbose = 0))
  expect_length(out, 0)
})

test_that("gl2geno output-file message names both real files at verbose >= 1", {
  # Expectation updated per approved finding F2: the message now prints the
  # two real output paths instead of the garbled '<outfile>.geno.lfmm.'.
  od <- file.path(tempdir(), "gl2geno_msg")
  dir.create(od, showWarnings = FALSE)
  out <- capture.output(gl2geno(testset2.gl, outfile = "msg", outpath = od,
                                verbose = 1))
  expect_true(any(grepl(file.path(od, "msg.geno"), out, fixed = TRUE)))
  expect_true(any(grepl(file.path(od, "msg.lfmm"), out, fixed = TRUE)))
})
