# Step 3, count from disk: report and filter functions read FBM-backed
# objects in blocks of loci (.gl_colwise, .gl_rowwise, .gl_n_na) instead of
# decoding the whole matrix. Results and messages must be identical to those
# of the gen-backed object, for any block size.

counts_pair <- function() {
  x <- gl.compliance.check(testset.gl[1:30, 1:60], verbose = 0)
  m <- as.matrix(x)
  m[, 5] <- NA                               # an all-NA locus
  m[, 6] <- 0L                               # a monomorphic locus
  m[pop(x) == levels(pop(x))[1], 7] <- NA    # all-NA in one population
  m[2, sample(ncol(m), 30)] <- NA            # a poorly scored individual
  x@gen <- new("genlight", m, ploidy = ploidy(x))@gen
  x <- gl.recalc.metrics(x, verbose = 0)
  list(x = x, f = suppressWarnings(gl.gen2fbm(x, verbose = 0)))
}

# The value and the console output of one call, with the object passed as
# 'x' in both storage types (history records the call)
run <- function(FUN, x, ...) {
  set.seed(1)                                # bootstraps draw loci
  out <- utils::capture.output(val <- suppressWarnings(FUN(x, ...)))
  if (methods::is(val, "genlight")) {
    val <- list(m = as.matrix(val), ind = indNames(val), loc = locNames(val),
                lm = val@other$loc.metrics)
  }
  if (inherits(val, "ggplot")) val <- NULL
  if (is.list(val) && !is.data.frame(val)) {
    val[vapply(val, inherits, logical(1), "ggplot")] <- NULL
  }
  list(val = val, out = out)
}

calls <- list(
  allele.freq = list(gl.allele.freq, simple = TRUE),
  allele.freq.pop = list(gl.allele.freq, by = "pop"),
  callrate.ind = list(gl.filter.callrate, method = "ind", threshold = 0.9,
                      plot.display = FALSE, verbose = 3),
  callrate.ind.rec = list(gl.filter.callrate, method = "ind",
                          threshold = 0.9, recursive = TRUE,
                          plot.display = FALSE, verbose = 3),
  report.callrate.loc = list(gl.report.callrate, method = "loc",
                             plot.display = FALSE, verbose = 3),
  report.callrate.ind = list(gl.report.callrate, method = "ind",
                             plot.display = FALSE, verbose = 3),
  monomorphs = list(gl.report.monomorphs, verbose = 3),
  maf = list(gl.report.maf, plot.display = FALSE, verbose = 3),
  rdepth = list(gl.report.rdepth, plot.display = FALSE, verbose = 3),
  reproducibility = list(gl.report.reproducibility, plot.display = FALSE,
                         verbose = 3),
  taglength = list(gl.report.taglength, plot.display = FALSE, verbose = 3),
  allna = list(gl.report.allna, verbose = 3),
  allna.pop = list(gl.report.allna, by.pop = TRUE, verbose = 3),
  heterozygosity = list(gl.filter.heterozygosity, t.upper = 0.3,
                        verbose = 3),
  basics = list(gl.report.basics, verbose = 2),
  tree.nj = list(gl.tree.nj, plot.display = FALSE, verbose = 3),
  het.pop = list(gl.report.heterozygosity, plot.display = FALSE,
                 verbose = 3),
  het.pop.inv = list(gl.report.heterozygosity, n.invariant = 50,
                     plot.display = FALSE, verbose = 3),
  het.pop.boot = list(gl.report.heterozygosity, nboots = 20,
                      CI.type = "perc", plot.display = FALSE, verbose = 3),
  het.ind = list(gl.report.heterozygosity, method = "ind",
                 plot.display = FALSE, verbose = 3),
  fst.pop = list(gl.fst.pop, verbose = 3),
  fst.pop.boot = list(gl.fst.pop, nboots = 5, verbose = 3),
  basic.stats = list(function(x, ...) {
    utils.basic.stats(gl.keep.pop(x, names(which(table(pop(x)) > 1)),
                                  verbose = 0), ...)
  }, rounded = FALSE),
  fstat = list(gl.report.fstat, plot.display = FALSE, verbose = 3),
  fstat.boot = list(gl.report.fstat, nboots = 10, CI.type = "perc",
                    plot.display = FALSE, verbose = 3),
  report.hwe = list(gl.report.hwe, sig_only = FALSE, min_sample_size = 2,
                    plot.out = FALSE, verbose = 3),
  filter.hwe = list(gl.filter.hwe, n.min = 2, verbose = 3),
  filter.hwe.all = list(gl.filter.hwe, subset = "all", verbose = 3),
  hwe.pop = list(gl.hwe.pop, plot.out = FALSE, HWformat = TRUE, verbose = 3),
  report.pa = list(gl.report.pa, loc.names = TRUE, matrix.pa = TRUE,
                   verbose = 3),
  report.pa.asym = list(gl.report.pa, test.asym = TRUE, test.asym.boot = 5,
                        verbose = 3),
  report.pa.1r = list(gl.report.pa, method = "one2rest", loc.names = TRUE,
                      verbose = 3),
  filter.pa = list(function(x, ...) {
    gl.filter.pa(x, popNames(x)[1], popNames(x)[2], ...)
  }, verbose = 3),
  filter.pa.inv = list(function(x, ...) {
    gl.filter.pa(x, popNames(x)[1], popNames(x)[2], ...)
  }, invers = TRUE, verbose = 3),
  diversity = list(gl.report.diversity, plot.display = FALSE, verbose = 3),
  diversity.1pop = list(function(x, ...) {
    gl.report.diversity(gl.keep.pop(x, popNames(x)[1], verbose = 0), ...)
  }, plot.display = FALSE, verbose = 3),
  allelerich = list(gl.report.allelerich, plot.display = FALSE, verbose = 3),
  allelerich.boot = list(gl.report.allelerich, nboots = 5, CI.type = "perc",
                         plot.display = FALSE, verbose = 3),
  test.het = list(gl.test.heterozygosity, nreps = 5, plot.out = FALSE,
                  verbose = 3),
  test.het.loc = list(gl.test.heterozygosity, nreps = 5, boot.method = "loc",
                      plot.out = FALSE, verbose = 3),
  shannon = list(gl.report.shannon, plot.display = FALSE, verbose = 3),
  shannon.beta = list(gl.report.shannon, level = "beta",
                      plot.display = FALSE, verbose = 3),
  shannon.gamma = list(gl.report.shannon, level = "gamma",
                       plot.display = FALSE, verbose = 3),
  fdsim = list(function(x, ...) {
    gl.fdsim(x, poppair = popNames(x)[1:2], ...)
  }, reps = 10, verbose = 3),
  filter.maf.pop = list(gl.filter.maf, by.pop = TRUE, ind.limit = 2,
                        pop.limit = 1, threshold = 0.1, verbose = 3),
  dist.ind = list(gl.dist.ind, plot.display = FALSE, verbose = 3),
  dist.ind.scaled = list(gl.dist.ind, scale = TRUE, type = "matrix",
                         plot.display = FALSE, verbose = 3),
  dist.ind.simple = list(gl.dist.ind, method = "Simple",
                         plot.display = FALSE, verbose = 3),
  dist.ind.absolute = list(gl.dist.ind, method = "Absolute",
                           plot.display = FALSE, verbose = 3),
  dist.ind.manhattan = list(gl.dist.ind, method = "Manhattan",
                            type = "matrix", plot.display = FALSE,
                            verbose = 3),
  dist.pop = list(gl.dist.pop, plot.display = FALSE, verbose = 3),
  dist.pop.nei = list(gl.dist.pop, method = "nei", plot.display = FALSE,
                      verbose = 3),
  replicates = list(gl.report.replicates, loc_threshold = 20,
                    perc_geno = 0.6, plot.out = FALSE, verbose = 3),
  replicates.king = list(gl.report.replicates.KING, loc_threshold = 20,
                         plot.out = FALSE, verbose = 3),
  contamination = list(gl.report.contamination, min.n = 2,
                       plot.display = FALSE, verbose = 3)
)

test_that("reports and filters give identical results for gen and FBM", {
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  old <- options(bigstatsr.block.sizeGB = getOption("bigstatsr.block.sizeGB"))
  on.exit(options(old), add = TRUE)
  set.seed(3)
  xs <- counts_pair()
  for (bs in c(1, 1e-6)) {                   # 1e-6: one locus per block
    options(bigstatsr.block.sizeGB = bs)
    for (nm in names(calls)) {
      cl <- calls[[nm]]
      g <- do.call(run, c(list(cl[[1]], xs$x), cl[-1]))
      f <- do.call(run, c(list(cl[[1]], xs$f), cl[-1]))
      expect_identical(f$val, g$val, label = paste(nm, bs))
      expect_identical(f$out, g$out, label = paste(nm, bs, "output"))
    }
  }
})

test_that("reports and filters never decode the whole FBM", {
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  set.seed(3)
  xs <- counts_pair()
  local_mocked_bindings(.fbm_as_matrix = function(x) {
    stop("the whole FBM was decoded")
  })
  for (nm in names(calls)) {
    cl <- calls[[nm]]
    # the function alone, without run(), whose comparison decodes the result
    set.seed(1)
    expect_no_error(utils::capture.output(suppressWarnings(
      do.call(cl[[1]], c(list(xs$f), cl[-1])))), message = nm)
  }
})

test_that(".gl_count_mean divides counts as rowMeans() does", {
  # pairs for which num / den differs from the long-double mean by one
  # unit in the last place (found by search)
  num <- c(3337, 3323, 14896, 0, 7, 5, NA, 0)
  den <- c(5390, 4590, 29237, 5000, 7, 10, 4, 0)
  ref <- vapply(seq_along(num), function(i) {
    if (is.na(num[i])) return(NA_real_)
    rowMeans(matrix(rep(c(TRUE, FALSE), c(num[i], den[i] - num[i])), 1),
             na.rm = TRUE)
  }, numeric(1))
  expect_identical(.gl_count_mean(num, den), ref)
})

test_that(".gl_pair_counts counts pairwise comparisons for any block size", {
  skip_if_not_installed("bigstatsr")
  skip_if_not_installed("Rcpp")
  old <- options(bigstatsr.block.sizeGB = getOption("bigstatsr.block.sizeGB"))
  on.exit(options(old), add = TRUE)
  set.seed(3)
  xs <- counts_pair()
  m <- as.matrix(xs$x)
  n <- nrow(m)
  ref <- list(L = integer(0), n1 = integer(0), n2 = integer(0))
  for (j in seq_len(n - 1)) for (i in (j + 1):n) {
    d <- abs(m[i, ] - m[j, ])
    ref$L <- c(ref$L, sum(!is.na(d)))
    ref$n1 <- c(ref$n1, sum(d == 1, na.rm = TRUE))
    ref$n2 <- c(ref$n2, sum(d == 2, na.rm = TRUE))
  }
  for (bs in c(1, 1e-6)) {
    options(bigstatsr.block.sizeGB = bs)
    expect_identical(.gl_pair_counts(xs$f), ref)
  }
  expect_error(.gl_pair_counts(xs$x), "file-backed")
  # one individual: no pairs
  expect_identical(.gl_pair_counts(xs$f[1, ]),
                   list(L = integer(0), n1 = integer(0), n2 = integer(0)))
})

test_that("exported utils.fbm.* helpers give the internal helpers' results", {
  skip_if_not_installed("bigstatsr")
  set.seed(3)
  xs <- counts_pair()
  r <- which(pop(xs$x) == levels(pop(xs$x))[1])
  for (x in xs) {
    cm <- function(m) colMeans(m, na.rm = TRUE)
    rs <- function(m) rowSums(is.na(m))
    expect_identical(utils.fbm.has(x), .has_fbm(x))
    expect_identical(utils.fbm.colwise(x, cm, rows = r),
                     .gl_colwise(x, cm, rows = r))
    expect_identical(utils.fbm.colblocks(x, cm), .gl_colblocks(x, cm))
    expect_identical(utils.fbm.rowwise(x, rs), .gl_rowwise(x, rs))
    expect_identical(utils.fbm.blocks(x), .fbm_col_blocks(nLoc(x), nInd(x)))
    expect_identical(utils.fbm.reader(x)(r, 3:9),
                     as.matrix(xs$x)[r, 3:9, drop = FALSE])
  }
  expect_false(utils.fbm.has(xs$x))
  expect_true(utils.fbm.has(xs$f))
  expect_identical(utils.fbm.pair.counts(xs$f), .gl_pair_counts(xs$f))
  expect_identical(utils.fbm.count.mean(3337, 5390),
                   .gl_count_mean(3337, 5390))
})
