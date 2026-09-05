# Characterization tests for gl2bpp — snapshot of CURRENT behaviour
# (review baseline, function-review campaign). Bugs are captured as-is and
# annotated; do not treat every expectation here as intended behaviour.

# Deterministic input used across tests (same pipeline as the roxygen example)
bpp_test_object <- function() {
  set.seed(42)
  test <- gl.filter.callrate(testset2.gl, threshold = 1, verbose = 0)
  test <- gl.filter.monomorphs(test, verbose = 0)
  gl.subsample.loc(test, n = 50, verbose = 0)
}

test_that("gl2bpp method 1 writes phylip-style locus blocks and an Imap", {
  od <- file.path(tempdir(), "bpp_m1")
  dir.create(od, showWarnings = FALSE)
  test <- bpp_test_object()

  ret <- gl2bpp(test, outpath = od, verbose = 0)
  expect_null(ret)

  bp <- readLines(file.path(od, "output_bpp.txt"))
  # 50 loci went in; gl.filter.overshoot inside gl2bpp removed 5, leaving 45
  # blocks of (1 header + 274 sequence lines)
  expect_length(bp, 45 * (nInd(test) + 1))
  expect_equal(trimws(bp[1]), "274 69")
  # observed 2026-09-05 (seed 42): first sequence line
  expect_equal(
    bp[2],
    paste0("SIM021-0-A/T^AA010915 ",
           "AAAAGGTTTTCATTGTACACCTCTGAATGAGACTCACCACGCACTCCGAGTGTCCTTAACCGTCCCTCT"))
  # [approved F6] lines are emitted with sep = "" -- no leading space on
  # sequence lines after the first in each block
  expect_false(startsWith(bp[3], " "))

  im <- readLines(file.path(od, "Imap.txt"))
  expect_length(im, nInd(test))
  expect_equal(trimws(im[1]), "AA010915 EmmacMDBForb")
})

test_that("gl2bpp method 2 writes a file of the same shape (seeded)", {
  od <- file.path(tempdir(), "bpp_m2")
  dir.create(od, showWarnings = FALSE)
  test <- bpp_test_object()
  set.seed(42)
  gl2bpp(test, method = 2, outfile = "m2.txt", outpath = od, verbose = 0)
  b2 <- readLines(file.path(od, "m2.txt"))
  expect_length(b2, 45 * (nInd(test) + 1))
  expect_equal(trimws(b2[1]), "274 69")
})

test_that("gl2bpp all-missing loci become all-N sequences for every individual", {
  od <- file.path(tempdir(), "bpp_n")
  dir.create(od, showWarnings = FALSE)
  # first 6 individuals of testset2.gl are all NA at locus 100049698-16-G/A
  a <- strsplit(as.character(testset2.gl@other$loc.metrics$AlleleID), "\\|")
  cloneid <- unlist(lapply(a, "[", 1))
  keep <- which(cloneid %in% c("100049687", "100049698", "100049728"))
  tsec <- testset2.gl[1:6, keep]
  gl2bpp(tsec, outfile = "n.txt", outpath = od, verbose = 0)
  bs <- readLines(file.path(od, "n.txt"))
  nblock <- grep("100049698-16-G/A", bs, value = TRUE)
  expect_length(nblock, 6)
  expect_true(all(grepl(" N+$", nblock)))
})

test_that("gl2bpp rejects SilicoDArT data", {
  expect_error(gl2bpp(testset2.gs, outpath = tempdir(), verbose = 0),
               "SilicoDArT")
})

test_that("gl2bpp rejects an invalid method with a fatal error", {
  # [approved F3] method is validated at entry; method = 3 previously ran
  # silently, wrote only the Imap, and reported Completed.
  od <- file.path(tempdir(), "bpp_m3")
  dir.create(od, showWarnings = FALSE)
  test <- bpp_test_object()
  expect_error(gl2bpp(test, method = 3, outfile = "m3.txt", outpath = od,
                      verbose = 0),
               "method must be")
  expect_false(file.exists(file.path(od, "m3.txt")))
})

test_that("gl2bpp runs when loc.metrics.flags is absent", {
  # [approved F4] absent flags are treated as monomorphs-not-confirmed
  # (gated warning); previously if(logical(0)) crashed with "argument is
  # of length zero".
  test <- bpp_test_object()
  test@other$loc.metrics.flags <- NULL
  expect_no_error(gl2bpp(test, outfile = "nf.txt", outpath = tempdir(),
                         verbose = 0))
})

test_that("gl2bpp output is identical for dartR-class and plain genlight input", {
  # [approved F1] the internal position sort now re-subsets loc.metrics
  # explicitly from the pre-sort object (DAT3), so plain genlight input no
  # longer pairs locus names with the wrong TrimmedSequence.
  od <- file.path(tempdir(), "bpp_cls")
  dir.create(od, showWarnings = FALSE)
  test <- bpp_test_object()
  gl2bpp(test, outfile = "dartR.txt", outpath = od, verbose = 0)
  pg <- test
  class(pg) <- "genlight"
  gl2bpp(pg, outfile = "genlight.txt", outpath = od, verbose = 0)
  d1 <- readLines(file.path(od, "dartR.txt"))
  d2 <- readLines(file.path(od, "genlight.txt"))
  expect_true(identical(d1, d2))
})

test_that("gl2bpp merge.secondaries merges clone blocks into one valid block", {
  # [approved F2] the merge block is rewritten: header recomputed from the
  # merged length, ALL superseded blocks deleted, non-overlapping segment
  # concatenation, labels indexed from the clone's first block. testset2.gl
  # carries no multi-locus clones, so a synthetic secondary pair is built.
  od <- file.path(tempdir(), "bpp_sec")
  dir.create(od, showWarnings = FALSE)
  sub <- testset2.gl[1:4, 1:3]
  tag <- paste(rep(c("A", "C", "G", "T"), 10), collapse = "")  # 40 bases
  lm <- sub@other$loc.metrics
  lm$AlleleID <- as.character(lm$AlleleID)
  lm$AlleleID[1] <- "900000001|F|0-5:C>T"
  lm$AlleleID[2] <- "900000001|F|0-20:G>A"
  lm$TrimmedSequence <- as.character(lm$TrimmedSequence)
  lm$TrimmedSequence[1] <- tag
  lm$TrimmedSequence[2] <- tag
  lm$SnpPosition[1] <- 5
  lm$SnpPosition[2] <- 20
  sub@other$loc.metrics <- lm
  sub@position[1] <- 5L
  sub@position[2] <- 20L
  sub@loc.all[1] <- "C/T"
  sub@loc.all[2] <- "G/A"

  gl2bpp(sub, outfile = "unmerged.txt", outpath = od, verbose = 0)
  gl2bpp(sub, outfile = "merged.txt", merge.secondaries = TRUE, outpath = od,
         verbose = 0)
  um <- readLines(file.path(od, "unmerged.txt"))
  mg <- readLines(file.path(od, "merged.txt"))

  # 3 loci in, 2 blocks out (the clone pair collapses to one merged block)
  expect_length(um, 3 * 5)
  expect_length(mg, 2 * 5)

  # the merged block replaces the clone's first (lowest-position) block and
  # carries a header recomputed from the merged sequence length
  sec <- grep("^sec-", mg)
  expect_length(sec, 4)
  expect_equal(mg[sec[1] - 1], "4 40")
  expect_true(all(startsWith(mg[sec], paste0("sec-", locNames(sub)[1], "^"))))

  # merged sequence = seq1[1..6] + seq2[7..21] + seq2[22..40]
  # (1-based SNP positions 6 and 21 for @position 5 and 20)
  blockseq <- function(lines, j) sub("^\\S+ ", "", lines[((j - 1) * 5 + 2):(j * 5)])
  s1 <- blockseq(um, 1)
  s2 <- blockseq(um, 2)
  expected <- paste0(substr(s1, 1, 6), substr(s2, 7, 21), substr(s2, 22, 40))
  expect_equal(sub("^\\S+ ", "", mg[sec]), expected)
  expect_true(all(nchar(expected) == 40))

  # the singleton locus block is untouched
  expect_equal(mg[6:10], um[11:15])
})

test_that("gl2bpp is silent at verbose = 0", {
  od <- file.path(tempdir(), "bpp_sil")
  dir.create(od, showWarnings = FALSE)
  test <- bpp_test_object()
  out <- capture.output(
    invisible(gl2bpp(test, outfile = "sil.txt", outpath = od, verbose = 0)))
  expect_length(out[out != "NULL"], 0)
})
