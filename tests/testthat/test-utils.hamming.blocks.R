# Characterization tests for utils.hamming.engine (utils.hamming.blocks.r)
# Baseline snapshotted before review (distance-wave review, dev at
# ddaed27). Engine verified sound; no flips.

test_that("dedup matches a brute-force worst-to-best scan", {
  skip_if_not_installed("Rcpp")
  eng <- utils.hamming.engine()
  set.seed(9)
  seqs <- lapply(1:40, function(i) as.raw(sample(1:4, 20, replace = TRUE)))
  brute <- rep(TRUE, 40)
  for (i in 40:1) {
    for (j in 40:1) {
      if (j > i && brute[j]) {
        mism <- sum(as.integer(seqs[[i]]) != as.integer(seqs[[j]]))
        if (mism <= 3) { brute[i] <- FALSE; break }
      }
    }
  }
  res <- eng$dedup(seqs, 3L)
  expect_identical(as.logical(res$keep), brute)
  expect_false(res$capped)
})

test_that("pairwise mismatch counts are exact", {
  skip_if_not_installed("Rcpp")
  eng <- utils.hamming.engine()
  set.seed(10)
  seqs <- lapply(1:12, function(i) as.raw(sample(1:4, 15, replace = TRUE)))
  pairs <- rbind(c(1L, 2L), c(3L, 10L), c(5L, 12L))
  pw <- eng$pairwise(seqs, pairs)
  manual <- apply(pairs, 1, function(p)
    sum(as.integer(seqs[[p[1]]]) != as.integer(seqs[[p[2]]])))
  expect_identical(as.integer(pw), as.integer(manual))
})
