# Characterization tests for gl.read.csv
# Baseline snapshotted before review (review-gl.read.csv, Phase A).
# Phase C (2026-09-05): [BUG] pins for the approved findings F1-F5 flipped
# to the fixed behaviour; each flipped expectation carries an
# "[approved F<n>]" comment. LOW findings (F6-F9) were deferred, so any
# behaviour they pin is unchanged.

# --- synthetic fixtures, written to tempdir ---------------------------------

write_fixture <- function(lines, name) {
  path <- file.path(tempdir(), name)
  writeLines(lines, path)
  path
}

num_basic <- write_fixture(c(
  "id,L1,L2,L3,L4,L5,L6",
  "i1,0,1,2,0,NA,1",
  "i2,1,1,0,0,0,2",
  "i3,2,0,1,1,0,0",
  "i4,0,2,NA,2,1,0",
  "i5,1,0,0,1,2,NA",
  "i6,2,2,1,0,1,0",
  "i7,0,0,2,NA,0,1",
  "i8,1,2,0,2,1,2"), "num_basic.csv")

# same file with an out-of-range genotype 5 at i3/L2
num_bad <- write_fixture(c(
  "id,L1,L2,L3,L4,L5,L6",
  "i1,0,1,2,0,NA,1",
  "i2,1,1,0,0,0,2",
  "i3,2,5,1,1,0,0",
  "i4,0,2,NA,2,1,0",
  "i5,1,0,0,1,2,NA",
  "i6,2,2,1,0,1,0",
  "i7,0,0,2,NA,0,1",
  "i8,1,2,0,2,1,2"), "num_bad.csv")

# character genotypes; per locus the designed ref allele is the most frequent:
# CL1 A/C (ref A), CL2 G/T (ref T), CL3 monomorphic A, CL4 C/T (ref C) with
# -/- missing, CL5 C/G (ref G); hets appear in both orders (A/C and C/A)
char_basic <- write_fixture(c(
  "id,CL1,CL2,CL3,CL4,CL5",
  "i1,A/A,T/T,A/A,C/C,G/G",
  "i2,A/C,T/T,A/A,-/-,G/G",
  "i3,C/A,G/T,A/A,C/T,G/G",
  "i4,A/A,T/T,A/A,T/T,G/C",
  "i5,C/C,T/G,A/A,C/C,G/G",
  "i6,A/A,G/G,A/A,-/-,G/G",
  "i7,A/C,T/T,A/A,C/C,C/G",
  "i8,A/A,T/T,A/A,T/C,G/G"), "char_basic.csv")

# char_basic with a third allele (G) at CL1/i8
char_tri <- write_fixture(c(
  "id,CL1,CL2,CL3,CL4,CL5",
  "i1,A/A,T/T,A/A,C/C,G/G",
  "i2,A/C,T/T,A/A,-/-,G/G",
  "i3,C/A,G/T,A/A,C/T,G/G",
  "i4,A/A,T/T,A/A,T/T,G/C",
  "i5,C/C,T/G,A/A,C/C,G/G",
  "i6,A/A,G/G,A/A,-/-,G/G",
  "i7,A/C,T/T,A/A,C/C,C/G",
  "i8,G/G,T/T,A/A,T/C,G/G"), "char_tri.csv")

small_num <- write_fixture(c(
  "id,L1,L2,L3",
  "i1,0,1,2",
  "i2,1,0,0",
  "i3,2,1,0"), "small_num.csv")

ind_meta_good <- write_fixture(c(
  "id,pop,sex",
  "i1,B,F", "i2,B,M", "i3,B,F", "i4,C,M",
  "i5,C,F", "i6,C,M", "i7,C,F", "i8,C,M"), "ind_meta_good.csv")

# identical content but id is the SECOND column
ind_meta_idsecond <- write_fixture(c(
  "pop,id,sex",
  "B,i1,F", "B,i2,M", "B,i3,F", "C,i4,M",
  "C,i5,F", "C,i6,M", "C,i7,F", "C,i8,M"), "ind_meta_idsecond.csv")

ind_meta_nopop <- write_fixture(c(
  "id,sex",
  "i1,F", "i2,M", "i3,F", "i4,M",
  "i5,F", "i6,M", "i7,F", "i8,M"), "ind_meta_nopop.csv")

loc_meta_good <- write_fixture(c(
  "AlleleID,SnpPosition,RepAvg",
  "L1,5,0.99", "L2,10,0.98", "L3,15,0.97",
  "L4,20,0.96", "L5,25,0.95", "L6,30,0.94"), "loc_meta_good.csv")

# --- numeric input -----------------------------------------------------------

test_that("numeric 0/1/2/NA csv reads faithfully", {
  o <- capture.output(g <- gl.read.csv(num_basic, verbose = 0))
  expect_length(o, 0)                       # verbose = 0 fully silent (VRB5)
  expect_s4_class(g, "dartR")
  expect_equal(nInd(g), 8)
  expect_equal(nLoc(g), 6)
  expect_equal(unique(ploidy(g)), 2)        # DAT1
  m <- as.matrix(g)
  raw <- as.matrix(read.csv(num_basic, row.names = 1))
  expect_equal(unname(m), unname(raw))      # every cell, NA included
  expect_equal(as.character(pop(g)), rep("A", 8))
  expect_equal(names(g@other$ind.metrics), c("id", "pop"))
  expect_equal(nrow(g@other$loc.metrics), 6)
  # genome-coordinate slots stay NULL (post-PR#330 convention)
  expect_null(g@position)
  expect_null(g@chromosome)
  expect_length(g@other$history, 1)
})

test_that("out-of-range numeric code 5 is a fatal error", {
  # [approved F2] validation now uses all-legal-codes logic; a file
  # containing 0,1,2 AND 5 errors, naming the offending value
  expect_error(gl.read.csv(num_bad, verbose = 0), "5")
})

# --- character (A/C/G/T pair) input -----------------------------------------

test_that("character genotypes convert with ref = most frequent allele", {
  g <- gl.read.csv(char_basic, verbose = 0)
  m <- as.matrix(g)
  exp <- rbind(c(0,0,0,0,0), c(1,0,0,NA,0), c(1,1,0,1,0), c(0,0,0,2,1),
               c(2,1,0,0,0), c(0,2,0,NA,0), c(1,0,0,0,1), c(0,0,0,1,0))
  # hand-verified against the raw file: hets in both orders (A/C, C/A) -> 1;
  # -/- -> NA; CL2 ref is T (most frequent) so G/G -> 2
  expect_equal(unname(m), exp)
  expect_equal(unique(ploidy(g)), 2)
})

test_that("real allele pairs are stored in loc.all", {
  # [approved F5] the character branch now stores the true ref/alt pair
  # per locus (ref = most frequent allele); the monomorphic CL3 carries
  # its single allele twice
  g <- gl.read.csv(char_basic, verbose = 0)
  expect_equal(g@loc.all, c("A/C", "T/G", "A/A", "C/T", "G/C"))
})

test_that("a third allele in character data is a fatal error", {
  expect_error(gl.read.csv(char_tri, verbose = 0))
})

test_that("transpose = TRUE recovers the same object", {
  d <- read.csv(char_basic, header = FALSE, stringsAsFactors = FALSE)
  tr <- file.path(tempdir(), "char_transposed.csv")
  write.table(t(d), tr, sep = ",", row.names = FALSE, col.names = FALSE,
              quote = FALSE)
  g1 <- gl.read.csv(char_basic, verbose = 0)
  g2 <- gl.read.csv(tr, transpose = TRUE, verbose = 0)
  expect_equal(as.matrix(g2), as.matrix(g1))
})

test_that("fewer than 5 loci reads correctly", {
  # [approved F3] the type-sniff window and the confirmation print are
  # clamped to the actual dimensions
  g <- gl.read.csv(small_num, verbose = 0)
  expect_equal(nInd(g), 3)
  expect_equal(nLoc(g), 3)
  m <- as.matrix(g)
  raw <- as.matrix(read.csv(small_num, row.names = 1))
  expect_equal(unname(m), unname(raw))
  # the verbose = 2 confirmation print no longer crashes either
  expect_no_error(capture.output(
    suppressWarnings(gl.read.csv(small_num, verbose = 2))))
})

# --- individual metadata -----------------------------------------------------

test_that("ind.metafile with id first joins and sets pop", {
  g <- gl.read.csv(num_basic, ind.metafile = ind_meta_good, verbose = 0)
  expect_equal(as.character(pop(g)), c("B","B","B","C","C","C","C","C"))
  expect_equal(as.character(g@other$ind.metrics$id), paste0("i", 1:8))
  expect_true("sex" %in% names(g@other$ind.metrics))
})

test_that("id column anywhere in the metafile is accepted", {
  # [approved F4] the order check now compares the id column by name;
  # a valid file with id second loads and sets pop correctly
  g <- gl.read.csv(num_basic, ind.metafile = ind_meta_idsecond,
                   verbose = 0)
  expect_equal(as.character(pop(g)), c("B","B","B","C","C","C","C","C"))
  expect_equal(as.character(g@other$ind.metrics$id), paste0("i", 1:8))
})

test_that("ind.metafile without pop column falls back to pop 'A'", {
  g <- suppressWarnings(
    gl.read.csv(num_basic, ind.metafile = ind_meta_nopop, verbose = 0))
  expect_equal(as.character(pop(g)), rep("A", 8))
  expect_true(all(c("id", "sex", "pop") %in% names(g@other$ind.metrics)))
})

# --- locus metadata ----------------------------------------------------------

test_that("a valid loc.metafile is validated and attached", {
  # [approved F1] the file's own AlleleID column is checked against the
  # locus names of the input data and the metrics are attached
  g <- gl.read.csv(num_basic, loc.metafile = loc_meta_good, verbose = 0)
  expect_true(all(c("AlleleID", "SnpPosition", "RepAvg") %in%
                    names(g@other$loc.metrics)))
  expect_equal(as.character(g@other$loc.metrics$AlleleID),
               paste0("L", 1:6))
  expect_equal(g@other$loc.metrics$RepAvg,
               c(0.99, 0.98, 0.97, 0.96, 0.95, 0.94))
})

test_that("a loc.metafile without AlleleID, or out of order, is a fatal error", {
  # [approved F1] the missing-column guard now stops (previously it
  # printed and fell through); order mismatches are reported by name
  loc_meta_noAID <- write_fixture(c(
    "SnpPosition,RepAvg",
    "5,0.99", "10,0.98", "15,0.97", "20,0.96", "25,0.95", "30,0.94"),
    "loc_meta_noAID.csv")
  expect_error(gl.read.csv(num_basic, loc.metafile = loc_meta_noAID,
                           verbose = 0), "AlleleID column absent")
  loc_meta_reorder <- write_fixture(c(
    "AlleleID,SnpPosition,RepAvg",
    "L2,10,0.98", "L1,5,0.99", "L3,15,0.97",
    "L4,20,0.96", "L5,25,0.95", "L6,30,0.94"), "loc_meta_reorder.csv")
  expect_error(gl.read.csv(num_basic, loc.metafile = loc_meta_reorder,
                           verbose = 0), "not in the same order")
  loc_meta_short <- write_fixture(c(
    "AlleleID,SnpPosition,RepAvg",
    "L1,5,0.99", "L2,10,0.98", "L3,15,0.97", "L4,20,0.96"),
    "loc_meta_short.csv")
  expect_error(gl.read.csv(num_basic, loc.metafile = loc_meta_short,
                           verbose = 0), "rows")
})

# --- fbm ---------------------------------------------------------------------

test_that("fbm = TRUE returns an FBM-backed object", {
  skip_if_not_installed("bigsnpr")
  g <- gl.read.csv(num_basic, fbm = TRUE, verbose = 0)
  expect_length(g@gen, 0)
  expect_false(is.null(g@fbm))
})
