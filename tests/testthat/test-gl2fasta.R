# Characterization tests for gl2fasta
# Baseline snapshotted before review (branch review-gl2fasta, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.
# Driven by a user report (Rachel, dartR group): a merged/co-analysed
# dataset yielded single-locus FASTA output with every method.

make_sub <- function() {
  sub <- gl.compliance.check(testset.gl[1:4, 1:10], verbose = 0)
  gl.filter.allna(sub, verbose = 0)
}
read_seqs <- function(f) {
  ln <- readLines(f)
  trimws(ln[!grepl("^>", ln)])
}

test_that("all four methods produce correct-length sequences", {
  skip_if_not_installed("seqinr")
  sub <- make_sub(); td <- tempdir()
  o <- capture.output({
    gl2fasta(sub, method = 1, outfile = "m1.fasta", outpath = td, verbose = 0)
    gl2fasta(sub, method = 2, outfile = "m2.fasta", outpath = td, verbose = 0)
    gl2fasta(sub, method = 3, outfile = "m3.fasta", outpath = td, verbose = 0)
    gl2fasta(sub, method = 4, outfile = "m4.fasta", outpath = td, verbose = 0)
  })
  expl <- sum(nchar(as.character(sub@other$loc.metrics$TrimmedSequence)))
  expect_equal(unique(nchar(read_seqs(file.path(td, "m1.fasta")))), expl)
  expect_equal(unique(nchar(read_seqs(file.path(td, "m2.fasta")))), expl)
  expect_equal(unique(nchar(read_seqs(file.path(td, "m3.fasta")))), nLoc(sub))
  expect_equal(unique(nchar(read_seqs(file.path(td, "m4.fasta")))), nLoc(sub))
  expect_equal(length(read_seqs(file.path(td, "m1.fasta"))), nInd(sub))
})

test_that("a gutting overshoot pre-filter is loud, not silent", {
  skip_if_not_installed("seqinr")
  # [approved diff F1] baseline: the internal gl.filter.overshoot call ran
  # at verbose 0; a dataset reduced to one locus produced single-locus
  # FASTA with no message at any verbosity (the user report).
  sub <- make_sub(); td <- tempdir()
  r <- sub
  r@other$loc.metrics$SnpPosition[2:nLoc(r)] <- 900
  expect_error(
    capture.output(gl2fasta(r, method = 3, outfile = "gut.fasta",
                            outpath = td, verbose = 0)),
    "overshoot")  # [approved diff F1]
})

test_that("non-numeric SnpPosition is fatal, not silently corrupting", {
  skip_if_not_installed("seqinr")
  # [approved diff F2] baseline: a factor SnpPosition ran substr with
  # factor LEVEL CODES, silently writing corrupt sequences (118/59/155
  # chars where 523 were expected on the fixture).
  sub <- make_sub(); td <- tempdir()
  # factor-coded numeric positions are recovered correctly (were corrupt)
  fct <- sub
  fct@other$loc.metrics$SnpPosition <- factor(fct@other$loc.metrics$SnpPosition)
  o <- capture.output(gl2fasta(fct, method = 1, outfile = "fct.fasta",
                               outpath = td, verbose = 0))
  o2 <- capture.output(gl2fasta(sub, method = 1, outfile = "ref.fasta",
                                outpath = td, verbose = 0))
  expect_identical(read_seqs(file.path(td, "fct.fasta")),
                   read_seqs(file.path(td, "ref.fasta")))  # [approved diff F2]
  # genuinely non-numeric positions are fatal
  bad <- sub
  bad@other$loc.metrics$SnpPosition <- as.character(bad@other$loc.metrics$SnpPosition)
  bad@other$loc.metrics$SnpPosition[1] <- "abc"
  expect_error(
    capture.output(gl2fasta(bad, method = 1, outfile = "bad.fasta",
                            outpath = td, verbose = 0)),
    "numeric")  # [approved diff F2]
})

test_that("method = 0 lists the options as documented", {
  sub <- make_sub()
  # [approved diff F4] baseline: the help path crashed on a cat():cat()
  # colon typo ("argument of length 0") before reaching its stop.
  expect_error(
    capture.output(gl2fasta(sub, method = 0, outpath = tempdir(),
                            verbose = 2)),
    "method")  # [approved diff F4]
})

test_that("return contract and flag robustness", {
  skip_if_not_installed("seqinr")
  sub <- make_sub(); td <- tempdir()
  o <- capture.output(rv <- gl2fasta(sub, method = 3, outfile = "rv.fasta",
                                     outpath = td, verbose = 0))
  expect_null(rv)
  nf <- sub
  nf@other$loc.metrics.flags <- NULL
  # [approved diff F5] baseline: the monomorphs check crashed on a
  # flag-less object ("argument is of length zero").
  expect_error(
    capture.output(gl2fasta(nf, method = 3, outfile = "nf.fasta",
                            outpath = td, verbose = 0)),
    NA)  # [approved diff F5]
})

test_that("an all-overshoot dataset stops cleanly without hijacking sinks", {
  skip_if_not_installed("seqinr")
  # [approved diff F1+F3] pre-fix this crashed mid-write inside an
  # unguarded sink, invalidating the session output stream.
  sub <- gl.compliance.check(testset.gl[1:4, 1:10], verbose = 0)
  sub <- gl.filter.allna(sub, verbose = 0)
  z <- sub
  z@other$loc.metrics$SnpPosition <- 900
  expect_error(
    capture.output(gl2fasta(z, method = 3, outfile = "z.fasta",
                            outpath = tempdir(), verbose = 0)),
    "overshoot")
  expect_equal(sink.number(), 0)
})
