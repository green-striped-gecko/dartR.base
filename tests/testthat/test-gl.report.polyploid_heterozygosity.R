# Characterization tests for gl.report.polyploid_heterozygosity, captured
# before the function-review changes and updated for the approved changes
# (report: function-review/reports/dartR.base/
# gl.report.polyploid_heterozygosity.md). Comments name the change that
# altered an assertion.

# Simulated autotetraploid dosage data (0-4 copies of the alternative allele),
# two populations of 20, 200 loci, 200 missing genotypes
sim_poly <- function(k = 4, n = 40, L = 200) {
  set.seed(42)
  q <- stats::runif(L, 0.05, 0.95)
  m <- sapply(q, function(qq) stats::rbinom(n, k, qq))
  m[sample(length(m), 200)] <- NA
  rownames(m) <- paste0("ind", seq_len(n))
  colnames(m) <- paste0("loc", seq_len(L))
  g <- new("genlight", m, ploidy = rep(k, n))
  pop(g) <- factor(rep(c("A", "B"), each = n / 2))
  list(gl = suppressWarnings(gl.compliance.check(g, verbose = 0)), m = m, k = k)
}

test_that("diploid pop output on testset.gl", {
  a <- gl.report.polyploid_heterozygosity(testset.gl, verbose = 0)
  expect_equal(dim(a), c(30L, 19L))
  expect_identical(names(a)[c(1, 8, 11, 14, 17)],
                   c("pop", "Ho", "He", "uHe", "FIS"))
  expect_equal(a$Ho[1:2], c(0.013560, 0.012813))
  expect_equal(a$He[1:2], c(0.013418, 0.009404))
  # addendum A2: all-NA loci are no longer counted as polymorphic
  # (previously 21 / 224 / 10 and 30 / 205 / 20)
  expect_equal(a$polyLoc[1:2], c(11, 10))
  expect_equal(a$monoLoc[1:2], c(234, 225))
  expect_equal(a$all_NALoc[1:2], c(10, 20))
  m1 <- as.matrix(seppop(testset.gl)[[1]])
  scored <- colSums(!is.na(m1)) > 0
  q1 <- colMeans(m1, na.rm = TRUE) / 2
  expect_equal(a$polyLoc[1], sum(scored & q1 > 0 & q1 < 1))
  # change 1: every other column now equals the diploid function (uHe and
  # FIS previously used the population's mean sample size)
  b <- gl.report.heterozygosity(testset.gl, verbose = 0)
  cols <- setdiff(names(a), c("polyLoc", "monoLoc"))
  expect_equal(a[, cols], b[, cols])
})

test_that("tetraploid pop output", {
  s <- sim_poly()
  a <- suppressWarnings(
    gl.report.polyploid_heterozygosity(s$gl, verbose = 0))
  # change 1: gametic Ho, dosage allele frequencies, per-locus N/(N-1),
  # checked against an independent computation
  k <- s$k
  for (i in 1:2) {
    mm <- s$m[as.integer(pop(s$gl)) == i, ]
    q <- colSums(mm, na.rm = TRUE) / (k * colSums(!is.na(mm)))
    He <- 2 * q * (1 - q)
    N <- k * colSums(!is.na(mm))
    Ho <- colMeans(mm * (k - mm) / choose(k, 2), na.rm = TRUE)
    uHe <- He * N / (N - 1)
    expect_equal(a$Ho[i], round(mean(Ho), 6))
    expect_equal(a$He[i], round(mean(He), 6))
    expect_equal(a$uHe[i], round(mean(uHe), 6))
    expect_equal(a$FIS[i], round(mean(1 - Ho / uHe), 6))
  }
})

test_that("tetraploid ind output", {
  s <- sim_poly()
  pdf(NULL)
  i <- suppressWarnings(
    gl.report.polyploid_heterozygosity(s$gl, method = "ind", verbose = 1))
  dev.off()
  expect_equal(nrow(i), 40L)
  # change 2: gametic Ho per individual; homozygotes are dosage 0 and k
  k <- s$k
  expect_equal(i$Ho, unname(rowMeans(s$m * (k - s$m) / choose(k, 2),
                                     na.rm = TRUE)))
  expect_equal(i$f.hom.ref, unname(rowMeans(s$m == 0, na.rm = TRUE)))
  expect_equal(i$f.hom.alt, unname(rowMeans(s$m == k, na.rm = TRUE)))
  # diploid individual output unchanged and equal to the diploid function
  pdf(NULL)
  d1 <- gl.report.polyploid_heterozygosity(testset.gl, method = "ind",
                                           verbose = 1)
  d2 <- gl.report.heterozygosity(testset.gl, method = "ind", verbose = 1)
  dev.off()
  expect_equal(d1, d2)
})

test_that("bootstrap intervals", {
  s <- sim_poly(k = 2)
  set.seed(1)
  o <- suppressWarnings(gl.report.polyploid_heterozygosity(
    s$gl, nboots = 100, CI.type = "perc", verbose = 0))
  # change 3: loci are resampled; intervals bracket their estimates
  expect_true(all(o$HeLCI <= o$He & o$He <= o$HeHCI))
  expect_true(all(o$HoLCI <= o$Ho & o$Ho <= o$HoHCI))
  s4 <- sim_poly()
  set.seed(1)
  o4 <- suppressWarnings(gl.report.polyploid_heterozygosity(
    s4$gl, nboots = 100, CI.type = "perc", verbose = 0))
  expect_true(all(o4$HeLCI <= o4$He & o4$He <= o4$HeHCI))
  expect_true(all(o4$HoLCI <= o4$Ho & o4$Ho <= o4$HoHCI))
})

test_that("former crashes now run", {
  pdf(NULL)
  on.exit(dev.off())
  # change 4
  expect_s3_class(gl.report.polyploid_heterozygosity(
    testset.gl, method = "ind", verbose = 0), "data.frame")
  expect_s3_class(suppressWarnings(gl.report.polyploid_heterozygosity(
    testset.gl, plot.display = FALSE, plot.file = "zz", verbose = 3)),
    "data.frame")
  expect_s3_class(suppressWarnings(gl.report.polyploid_heterozygosity(
    testset.gl, method = "ind", subsample.pop = TRUE, verbose = 2)),
    "data.frame")
  expect_type(gl.report.polyploid_heterozygosity(
    testset.gl, subsample.pop = TRUE, verbose = 3), "list")
  set.seed(1)
  expect_s3_class(suppressWarnings(gl.report.polyploid_heterozygosity(
    testset.gl, nboots = 20, CI.type = "perc", verbose = 0)), "data.frame")
  # change 5: the error carries its message
  expect_error(gl.report.polyploid_heterozygosity(
    testset.gl, error.bar = "CI", verbose = 0), "nboots")
})

test_that("silent at verbose 0", {
  # change 5
  expect_silent(gl.report.polyploid_heterozygosity(testset.gl,
                                                   method = "bogus",
                                                   verbose = 0))
})

test_that("tetraploid subsample Ho is gametic", {
  # addendum A1: a subsample of every individual reproduces gametic Ho
  set.seed(42)
  k <- 4
  m <- sapply(stats::runif(100, 0.1, 0.9),
              function(qq) stats::rbinom(20, k, qq))
  g <- new("genlight", m, ploidy = rep(k, 20))
  pop(g) <- factor(rep("A", 20))
  r <- utils.subsample.pop(g, n.limit = 10, subsamples = 20)
  expect_equal(r$res.mean, mean(colMeans(m * (k - m) / choose(k, 2))))
})

test_that("input untouched, subsample return shape", {
  x0 <- testset.gl
  r <- gl.report.polyploid_heterozygosity(x0, subsample.pop = TRUE,
                                          verbose = 0)
  expect_identical(x0, testset.gl)
  # change 6: named list
  expect_named(r, c("subsample", "results"))
  expect_s3_class(r$results, "data.frame")
})
