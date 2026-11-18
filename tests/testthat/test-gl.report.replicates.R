# Characterization tests for gl.report.replicates
# Baseline snapshotted before review (review-gl.report.replicates).
# Assertions marked [approved diff] were flipped in Phase C.

skip_if_not_installed("Rcpp")
skip_if_not_installed("RcppParallel")

# platypus subset with two exact duplicate individuals injected
make_dup <- function() {
  base <- platypus.gl[1:30, 1:500]
  bm <- as.matrix(base)
  dm <- rbind(bm, bm[1:2, ])
  rownames(dm) <- make.unique(c(rownames(bm), rownames(bm)[1:2]))
  dup <- new("genlight", gen = dm, ploidy = 2)
  indNames(dup) <- rownames(dm)
  locNames(dup) <- colnames(bm)
  pop(dup) <- factor(rep("A", nrow(dm)))
  suppressWarnings(gl.compliance.check(dup, verbose = 0))
}
dup <- make_dup()

test_that("pair statistics verified against direct computation", {
  pdf(NULL); on.exit(dev.off())
  o <- capture.output(res <- gl.report.replicates(dup, loc_threshold = 100,
        perc_geno = 0.9, plot.out = FALSE, verbose = 0))
  expect_length(o, 0)
  expect_type(res, "list")
  expect_named(res, c("table.rep", "ind.list.drop", "ind.list.rep"))
  tr <- res$table.rep
  # ind1 is the individual earlier in the object (post-F1 orientation)
  row <- tr[tr$ind1 == "T3" & tr$ind2 == "T5", ]
  expect_equal(row$nloc, 468)
  expect_equal(round(row$perc, 5), 0.99359)
  # independent recomputation of that pair
  m <- as.matrix(dup)
  bothloc <- !is.na(m["T5", ]) & !is.na(m["T3", ])
  expect_equal(sum(bothloc), 468)
  expect_equal(round(sum(m["T5", bothloc] == m["T3", bothloc]) / 468, 5),
               round(row$perc, 5))
})

test_that("duplicate pairs and tied missingness", {
  pdf(NULL); on.exit(dev.off())
  invisible(capture.output(res <- gl.report.replicates(dup,
        loc_threshold = 100, perc_geno = 0.9, plot.out = FALSE,
        verbose = 0)))
  # [approved diff F1] baseline: every pair appeared twice in table.rep
  # (both orderings; 3 real pairs -> 6 rows), and for tied missingness
  # the drop rule selected the OPPOSITE member in each ordering, putting
  # BOTH members of an exact-duplicate pair in ind.list.drop.
  expect_equal(nrow(res$table.rep), 3)
  expect_setequal(res$ind.list.drop, c("T27.1", "T35.1", "T3"))
  # replicate partners listed for exactly the six involved individuals
  expect_setequal(names(res$ind.list.rep),
                  c("T27", "T27.1", "T35", "T35.1", "T3", "T5"))
})

test_that("no-pairs path", {
  pdf(NULL); on.exit(dev.off())
  base <- platypus.gl[1:30, 1:500]
  o <- capture.output(res <- gl.report.replicates(base, loc_threshold = 100,
        perc_geno = 0.999, plot.out = FALSE, verbose = 0))
  # [approved diff F2] baseline: returned a bare character message.
  expect_type(res, "list")
  expect_named(res, c("table.rep", "ind.list.drop", "ind.list.rep"))
  expect_equal(nrow(res$table.rep), 0)
  expect_length(res$ind.list.drop, 0)
  expect_length(o, 0)
})

test_that("return visibility and input integrity", {
  pdf(NULL); on.exit(dev.off())
  xcopy <- dup
  invisible(capture.output(v <- withVisible(gl.report.replicates(dup,
        loc_threshold = 100, perc_geno = 0.9, plot.out = FALSE,
        verbose = 0))))
  expect_false(v$visible)     # [approved diff F4] was returned visibly
  expect_identical(xcopy, dup)
})

test_that("plot rendering at verbose 0", {
  # BASELINE: no verbose==0 gate on plot.out; the plot renders at v0.
  t0 <- tempfile(fileext = ".pdf"); pdf(t0); dev.off()
  empty <- file.info(t0)$size
  tf <- tempfile(fileext = ".pdf"); pdf(tf)
  invisible(capture.output(r <- gl.report.replicates(dup,
        loc_threshold = 100, perc_geno = 0.9, verbose = 0)))
  dev.off()
  # [approved diff F3] baseline: a page was drawn at verbose 0.
  expect_lt(file.info(tf)$size, empty + 500)
  unlink(c(t0, tf))
})
