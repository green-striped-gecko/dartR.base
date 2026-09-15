# Characterization tests for gl.test.heterozygosity
# Baseline snapshotted before review (review-gl.test.heterozygosity) at
# 98d3a98. Assertions marked [approved diff, change n] were flipped in
# Phase C to reflect the approved behaviour changes (report:
# function-review/reports/dartR.base/gl.test.heterozygosity.md). Bootstraps are seeded; nreps is kept
# small for speed, so numerical pins are seed-specific.

indep_uHe <- function(x) {
  vapply(seppop(x), function(p) {
    m <- as.matrix(p)
    n <- colSums(!is.na(m))
    pf <- colMeans(m, na.rm = TRUE) / 2
    mean((2 * n / (2 * n - 1)) * 2 * pf * (1 - pf), na.rm = TRUE)
  }, numeric(1))
}

run_quiet <- function(...) {
  out <- capture.output(res <- gl.test.heterozygosity(..., verbose = 0))
  list(res = res, out = out)
}

test_that("platypus.gl, ind bootstrap: structure and observed differences", {
  set.seed(1)
  r <- run_quiet(platypus.gl, nreps = 50, plot.out = FALSE)
  res <- r$res
  expect_s3_class(res, "data.frame")
  expect_equal(dim(res), c(3, 8))
  expect_equal(colnames(res), c("pop1", "pop2", "diff", "CIlower", "CIupper",
                                "significance", "pval", "pval.adj"))
  expect_equal(res$pop1, c("SEVERN_ABOVE", "SEVERN_ABOVE", "SEVERN_BELOW"))
  expect_equal(res$pop2, c("SEVERN_BELOW", "TENTERFIELD", "TENTERFIELD"))
  # observed differences equal an independent per-locus uHe computation
  u <- indep_uHe(platypus.gl)
  expect_equal(unname(res$diff), unname(c(u[1] - u[2], u[1] - u[3], u[2] - u[3])),
               tolerance = 1e-10)
  # seed-specific pins
  expect_equal(res$CIlower, c(-0.0013850, -0.0146799, -0.0178369), tolerance = 1e-5)
  expect_equal(res$CIupper, c(0.00905783, -0.00619746, -0.00919149), tolerance = 1e-5)
  expect_equal(res$pval, c(0.313725, 0.0392157, 0.0392157), tolerance = 1e-5)
  expect_equal(res$pval.adj, p.adjust(res$pval, "BH"), tolerance = 1e-5)
  # [approved diff, change 1] labels now use alpha/2 tails; at nreps = 50 the
  # p value floor is 2/51 = 0.039 while all replicates of pairs 2-3 are < 0
  expect_equal(res$significance, c("non-sig @0.05", "sig @0.01", "sig @0.01"))
  # [approved diff, change 5] nothing printed at verbose = 0
  expect_length(r$out, 0)
  # read-only
  expect_identical(platypus.gl, dartR.data::platypus.gl)
})

test_that("F1: labels, p values and CI agree at nreps = 1000", {
  big <- names(sort(table(pop(testset.gl)), decreasing = TRUE))[1:8]
  x <- gl.keep.pop(testset.gl, pop.list = big, verbose = 0)
  set.seed(7)
  res <- run_quiet(x, nreps = 1000, plot.out = FALSE)$res
  sig <- grepl("^sig", res$significance)
  # [approved diff, change 1] was 2 pairs labelled sig @0.05 with p 0.07-0.09
  expect_equal(sum(sig & res$pval > 0.05), 0)
  expect_equal(sum(!sig & res$pval <= 0.05), 0)
  expect_equal(sum(res$significance == "sig @0.01" & res$pval > 0.01), 0)
  expect_equal(sum(sig & res$CIlower < 0 & res$CIupper > 0), 0)
  expect_equal(sum(res$significance == "sig @0.01"), 10)
})

test_that("PLT3: results independent of plotting; plot files written", {
  set.seed(1)
  r0 <- run_quiet(platypus.gl, nreps = 50, plot.out = FALSE)$res
  td <- tempfile(); dir.create(td)
  set.seed(1)
  pdf(NULL)
  r1 <- run_quiet(platypus.gl, nreps = 50, plot.out = TRUE,
                  plot.file = "het", plot.dir = td)$res
  dev.off()
  expect_identical(r0, r1)
  # one RDS per page of plots plus the table RDS (documented, change 7)
  expect_setequal(list.files(td), c("het_1_to_3.RDS", "table_het.RDS"))
  expect_s3_class(readRDS(file.path(td, "het_1_to_3.RDS")), "patchwork")
  expect_identical(readRDS(file.path(td, "table_het.RDS")), r0)
})

test_that("boot.method = 'loc' resamples loci", {
  set.seed(1)
  res <- run_quiet(platypus.gl, nreps = 50, boot.method = "loc",
                   plot.out = FALSE)$res
  expect_equal(res$diff, c(0.002007828, -0.008589002, -0.010596830),
               tolerance = 1e-6)
  expect_equal(res$CIlower, c(-0.0126488, -0.0179299, -0.0206603), tolerance = 1e-5)
  expect_equal(res$significance, rep("non-sig @0.05", 3))
})

test_that("paging: 5 populations, 10 pairs, max_plots 6", {
  big <- names(sort(table(pop(testset.gl)), decreasing = TRUE))[1:5]
  x <- gl.keep.pop(testset.gl, pop.list = big, verbose = 0)
  set.seed(1)
  pdf(NULL)
  res <- run_quiet(x, nreps = 10, plot.out = TRUE)$res
  dev.off()
  expect_equal(nrow(res), 10)
})

test_that("FBM and SilicoDArT inputs", {
  xf <- gl.gen2fbm(platypus.gl, verbose = 0)
  set.seed(1)
  rf <- run_quiet(xf, nreps = 20, plot.out = FALSE)$res
  set.seed(1)
  rd <- run_quiet(platypus.gl, nreps = 20, plot.out = FALSE)$res
  expect_identical(rf, rd)
  # [approved diff, change 2] SilicoDArT is rejected
  xs <- gl.keep.pop(testset.gs, pop.list = popNames(testset.gs)[1:3], verbose = 0)
  expect_error(run_quiet(xs, nreps = 20, plot.out = FALSE),
               "found SilicoDArT expecting SNP")
})

test_that("F3/F4: fewer than two populations error early; missing flags tolerated", {
  # [approved diff, change 3] was "subscript out of bounds" after the bootstrap
  x1 <- gl.keep.pop(platypus.gl, pop.list = popNames(platypus.gl)[1], verbose = 0)
  expect_error(run_quiet(x1, nreps = 10, plot.out = FALSE),
               "At least two populations")
  x0 <- platypus.gl; pop(x0) <- NULL
  expect_error(run_quiet(x0, nreps = 10, plot.out = FALSE),
               "At least two populations")
  # [approved diff, change 4] was "argument is of length zero"
  g <- new("genlight", as.matrix(platypus.gl)[1:30, 1:50],
           pop = pop(platypus.gl)[1:30])
  set.seed(1)
  rg <- run_quiet(g, nreps = 10, plot.out = FALSE)$res
  expect_equal(nrow(rg), 3)
  u <- indep_uHe(g)
  expect_equal(unname(rg$diff), unname(c(u[1] - u[2], u[1] - u[3], u[2] - u[3])),
               tolerance = 1e-10)
})

test_that("F5/F6: warnings gated by verbose; alpha order", {
  set.seed(1)
  r <- run_quiet(platypus.gl, nreps = 10, alpha1 = 2, boot.method = "x",
                 plot.out = FALSE)
  # [approved diff, change 5] warnings silent at verbose = 0, shown at 2
  expect_length(r$out, 0)
  set.seed(1)
  out <- capture.output(gl.test.heterozygosity(platypus.gl, nreps = 10,
                                               alpha1 = 2, boot.method = "x",
                                               plot.out = FALSE, verbose = 2))
  expect_true(any(grepl("First alpha value", out)))
  expect_true(any(grepl("boot.method must be", out)))
  expect_equal(out[1], "Starting gl.test.heterozygosity ")
  # alpha1 > alpha2: table labels use the looser level for @0.05 ...
  set.seed(1)
  res <- run_quiet(platypus.gl, nreps = 50, alpha1 = 0.01, alpha2 = 0.05,
                   plot.out = FALSE)$res
  expect_equal(res$significance, c("non-sig @0.05", "sig @0.01", "sig @0.01"))
  # ... and [approved diff, change 6] the legend labels follow the swap
  td <- tempfile(); dir.create(td)
  set.seed(1)
  pdf(NULL)
  run_quiet(platypus.gl, nreps = 20, alpha1 = 0.01, alpha2 = 0.05,
            plot.out = TRUE, plot.file = "het", plot.dir = td)
  dev.off()
  p <- readRDS(file.path(td, "het_1_to_3.RDS"))
  lab <- ggplot2::ggplot_build(p[[3]])$plot$scales$get_scales("colour")$labels
  expect_equal(lab, c("Sig.  0.05", "Sig.  0.01", "Observed", "Zero value"))
})
