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
