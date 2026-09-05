# Characterization tests for gl2treemix
# Baseline snapshotted pre-review (dev at ed99203). Assertions tagged
# [bug baseline] capture current defective behaviour on purpose; flip them
# only against an approved finding in the Phase A report.

tm_sub <- function() testset2.gl[1:8, 1:15]

test_that("gzipped treemix file has header plus one aligned row per locus", {
  td <- tempdir()
  f <- file.path(td, "tm_base.gz")
  o <- capture.output(r <- gl2treemix(tm_sub(), outfile = "tm_base.gz",
                                      outpath = td, verbose = 0))
  expect_null(r)
  con <- gzfile(f)
  ln <- readLines(con)
  close(con)
  expect_length(ln, 1 + nLoc(tm_sub()))
  hdr <- strsplit(trimws(ln[1]), " +")[[1]]
  # header lists the populations present, in sorted factor-level order
  expect_setequal(hdr, as.character(unique(pop(tm_sub()))))
  # every cell must equal the manual diploid allele count refN,altN for
  # its (locus, population); verifies row order and column order jointly
  gm <- as.matrix(tm_sub())
  pops <- as.character(pop(tm_sub()))
  for (j in seq_len(nLoc(tm_sub()))) {
    cells <- strsplit(trimws(ln[j + 1]), " +")[[1]]
    expect_length(cells, length(hdr))
    for (pi in seq_along(hdr)) {
      v <- gm[pops == hdr[pi], j]
      alt <- sum(v, na.rm = TRUE)
      ref <- 2 * sum(!is.na(v)) - alt
      expect_equal(cells[pi], paste0(ref, ",", alt),
                   label = sprintf("locus %d pop %s", j, hdr[pi]))
    }
  }
  # a population scored all-NA at a locus yields the treemix missing code 0,0
  expect_true(any(unlist(strsplit(trimws(ln[-1]), " +")) == "0,0"))
})

test_that("SilicoDArT is rejected explicitly", {
  # F1 fix applied: accept = 'SNP' (DAT7) makes ploidy-1 data a fatal
  # error instead of silently doubling allele copies
  expect_error(
    capture.output(gl2treemix(testset2.gs[1:6, 1:5], outfile = "tm_gs.gz",
                              outpath = tempdir(), verbose = 0)),
    "SilicoDArT")
})

test_that("no connection is left open after writing", {
  td <- tempdir()
  before <- nrow(showConnections())
  capture.output(gl2treemix(tm_sub(), outfile = "tm_con.gz", outpath = td,
                            verbose = 0))
  expect_equal(nrow(showConnections()), before)
  expect_equal(sink.number(), 0)
})

test_that("NULL return is invisible", {
  # F4 fix applied: invisible(NULL) rather than return(NULL)
  o <- capture.output(
    vis <- withVisible(gl2treemix(tm_sub(), outfile = "tm_vis.gz",
                                  outpath = tempdir(), verbose = 0)))
  expect_false(vis$visible)
})
