# Characterization tests for gl.report.diversity
# Baseline snapshotted before review (branch review-gl.report.diversity,
# dev at ddaed27). Assertions marked [approved diff] were flipped in
# Phase C.

make_g7 <- function() {
  set.seed(7)
  m <- matrix(sample(0:2, 12 * 20, replace = TRUE), nrow = 12)
  m[1:6, 3] <- NA   # locus 3 all-NA in pop P1
  g7 <- new("genlight", gen = m, ploidy = 2)
  indNames(g7) <- paste0("i", 1:12); locNames(g7) <- paste0("L", 1:20)
  pop(g7) <- factor(rep(c("P1", "P2"), each = 6))
  suppressWarnings(gl.compliance.check(g7, verbose = 0))
}

shan_vec <- function(gm) {
  A <- 2 * colSums(gm == 0, na.rm = TRUE) + colSums(gm == 1, na.rm = TRUE)
  B <- 2 * colSums(gm == 2, na.rm = TRUE) + colSums(gm == 1, na.rm = TRUE)
  sapply(seq_along(A), function(i) {
    v <- c(A[i], B[i]); v <- v[v > 0]
    if (!length(v)) return(NA)
    pr <- v / sum(v); -sum(pr * log(pr))
  })
}

test_that("q=0 and q=2 alpha match an independent recomputation", {
  xf <- gl.filter.allna(testset.gl, verbose = 0)
  o <- capture.output(r <- gl.report.diversity(testset.gl,
        plot.display = FALSE, verbose = 0))
  gm <- as.matrix(xf[pop(xf) == popNames(xf)[1], ])
  cm <- colMeans(gm, na.rm = TRUE)
  expect_equal(unname(r$zero_H_alpha[1]), mean((cm %% 2) > 0, na.rm = TRUE))
  p <- cm / 2; p <- p[!is.na(p)]
  expect_equal(unname(r$two_H_alpha[1]), mean(2 * p * (1 - p)))
})

test_that("verbose 0 is silent", {
  # [approved diff F3] baseline: the kable tables printed unconditionally
  # (356 lines at verbose 0 on testset).
  o <- capture.output(r <- gl.report.diversity(testset.gl,
        plot.display = FALSE, verbose = 0))
  expect_equal(length(o), 0)  # [approved diff F3]
})

test_that("Shannon alpha excludes per-pop all-missing loci", {
  # [approved diff F1] baseline: all-NA loci entered the mean as zero
  # diversity, deflating one_H_alpha.
  g7 <- make_g7()
  o <- capture.output(r <- gl.report.diversity(g7, plot.display = FALSE,
        verbose = 0))
  s1 <- shan_vec(as.matrix(g7[pop(g7) == "P1", ]))
  expect_equal(unname(r$one_H_alpha[1]),
               mean(s1, na.rm = TRUE))  # [approved diff F1]
})

test_that("Shannon beta matches a pair-pooled recomputation", {
  # [approved diff F1] baseline: the full-length dummys vector misaligned
  # with the beta indexing (logical recycling) whenever a population
  # carried an all-NA locus.
  g7 <- make_g7()
  o <- capture.output(r <- gl.report.diversity(g7, plot.display = FALSE,
        verbose = 0))
  gm1 <- as.matrix(g7[pop(g7) == "P1", ])
  gm2 <- as.matrix(g7[pop(g7) == "P2", ])
  s1 <- shan_vec(gm1); s2 <- shan_vec(gm2)
  pall <- colMeans(rbind(gm1, gm2), na.rm = TRUE) / 2
  Hall <- -(pall * log(pall) + (1 - pall) * log(1 - pall))
  Hall[!is.finite(Hall)] <- 0
  idx <- which(!is.na(s1) & !is.na(s2))
  manual <- mean(Hall[idx] - (s1[idx] + s2[idx]) / 2)
  expect_equal(unname(r$one_H_beta[2, 1]), manual)  # [approved diff F1+F2]
})

test_that("pairwise betas depend only on the pair", {
  # [approved diff F2] baseline: one_H_beta and two_H_beta pooled ALL
  # populations (0.0396 vs 0.0103 and 0.0211 vs 0.0108 for the same pair
  # on testset); zero_H_beta was already pair-based.
  pn <- popNames(testset.gl)
  pair <- gl.keep.pop(testset.gl, pop.list = pn[1:2], verbose = 0)
  o1 <- capture.output(rfull <- gl.report.diversity(testset.gl,
        plot.display = FALSE, verbose = 0))
  o2 <- capture.output(rpair <- gl.report.diversity(pair,
        plot.display = FALSE, verbose = 0))
  expect_equal(rfull$zero_H_beta[2, 1], rpair$zero_H_beta[2, 1])
  expect_equal(rfull$one_H_beta[2, 1], rpair$one_H_beta[2, 1])  # [approved diff F2]
  expect_equal(rfull$two_H_beta[2, 1], rpair$two_H_beta[2, 1])  # [approved diff F2]
})

test_that("plot.file without plot.display works", {
  # [approved diff F4] baseline: "object 'p3' not found".
  expect_error(
    gl.report.diversity(testset.gl, plot.display = FALSE,
        plot.file = "divtest", plot.dir = tempdir(), table = "N",
        verbose = 0),
    NA)  # [approved diff F4]
})

test_that("SilicoDArT input is rejected", {
  # [approved diff F5] baseline: silico data ran silently through diploid
  # 0/1/2 arithmetic, producing meaningless indices.
  expect_error(
    capture.output(r <- gl.report.diversity(testset.gs,
        plot.display = FALSE, verbose = 0)),
    "SNP")  # [approved diff F5]
})

test_that("an invalid table argument fails informatively", {
  # [approved diff F3] baseline: a typo silently suppressed all tables.
  pn <- popNames(testset.gl)
  pair <- gl.keep.pop(testset.gl, pop.list = pn[1:2], verbose = 0)
  expect_error(
    capture.output(gl.report.diversity(pair, plot.display = FALSE,
        table = "HDX", verbose = 2)),
    "table")  # [approved diff F3]
})

test_that("plotting at verbose 0 is silent", {
  # [approved diff F6] baseline: the gl.colors banner leaked.
  pdf(NULL); on.exit(dev.off())
  o <- capture.output(r <- gl.report.diversity(testset.gl,
        plot.display = TRUE, table = "N", verbose = 0))
  expect_equal(length(o), 0)  # [approved diff F6]
})
