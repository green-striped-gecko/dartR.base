# Characterization tests for gl.randomize.snps, captured before the
# function-review changes and updated for the approved changes (report:
# function-review/reports/dartR.base/gl.randomize.snps.md). Comments name the
# change that altered an assertion.

x <- gl.filter.monomorphs(testset.gl, verbose = 0)

test_that("SNP recoding swaps 0 and 2 in half the loci", {
  set.seed(1)
  y <- gl.randomize.snps(x, plot.display = FALSE, verbose = 0)
  m0 <- as.matrix(x)
  m1 <- as.matrix(y)
  changed <- which(colSums(m0 != m1, na.rm = TRUE) > 0)

  expect_equal(dim(m1), dim(m0))
  expect_length(changed, floor(nLoc(x) / 2))
  expect_equal(unique(ploidy(y)), 2L)
  expect_identical(is.na(m1), is.na(m0))
  expect_identical(m1 == 1, m0 == 1)
  expect_equal(unname(glMean(y)[changed]), unname(1 - glMean(x)[changed]))
  expect_identical(locNames(y), locNames(x))
  expect_identical(indNames(y), indNames(x))
  expect_identical(pop(y), pop(x))
  expect_identical(nrow(y@other$loc.metrics), nLoc(x))
  expect_identical(nrow(y@other$ind.metrics), nInd(x))
  # change 3: allele labels reversed for recoded loci only
  flip <- function(a) vapply(strsplit(a, "/", fixed = TRUE),
                             function(z) paste(rev(z), collapse = "/"), "")
  expect_identical(y@loc.all[changed], flip(x@loc.all[changed]))
  expect_identical(y@loc.all[-changed], x@loc.all[-changed])
})

test_that("flags reset and history appended", {
  y <- gl.randomize.snps(x, plot.display = FALSE, verbose = 0)
  expect_false(any(unlist(y@other$loc.metrics.flags)))
  expect_length(y@other$history, length(x@other$history) + 1)
  expect_true(validObject(y))
})

test_that("plotting does not change the returned genotypes", {
  set.seed(1)
  y1 <- gl.randomize.snps(x, plot.display = FALSE, verbose = 0)
  set.seed(1)
  pdf(NULL)
  y2 <- gl.randomize.snps(x, plot.display = TRUE, verbose = 1)
  dev.off()
  expect_identical(as.matrix(y2), as.matrix(y1))
})

test_that("SilicoDArT input is refused", {
  # change 1: previously wrote invalid 2 codes into a ploidy-1 object
  expect_error(gl.randomize.snps(testset.gs, verbose = 0))
})

test_that("FBM-backed input", {
  skip_if_not_installed("bigstatsr")
  xf <- gl.gen2fbm(x, verbose = 0)
  set.seed(1)
  zf <- gl.randomize.snps(xf, plot.display = FALSE, verbose = 0)
  # change 2: the FBM copy is recoded, matching the dense result for the same
  # seed, and the input object is untouched
  set.seed(1)
  y <- gl.randomize.snps(x, plot.display = FALSE, verbose = 0)
  expect_equal(unname(as.matrix(zf)), unname(as.matrix(y)))
  expect_identical(zf@loc.all, y@loc.all)
  expect_equal(unname(as.matrix(xf)), unname(as.matrix(x)))
})

test_that("plot.colors and plot.theme reach the smear plots", {
  # change 6: previously accepted and ignored
  seen <- list()
  local_mocked_bindings(gl.smearplot = function(x, ..., plot.theme = NULL,
                                                plot.colors = NULL) {
    seen[[length(seen) + 1]] <<- list(theme = plot.theme, cols = plot.colors)
    ggplot2::ggplot()
  })
  cols <- c("black", "grey", "white", "pink")
  pdf(NULL)
  gl.randomize.snps(x, plot.colors = cols, plot.theme = ggplot2::theme_bw(),
                    verbose = 1)
  dev.off()
  expect_length(seen, 2)
  for (s in seen) {
    expect_identical(s$cols, cols)
    expect_s3_class(s$theme, "theme")
  }
})
