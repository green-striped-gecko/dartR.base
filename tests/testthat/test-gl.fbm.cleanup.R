# Slow file: skipped on CRAN to keep the check time under 10 minutes; runs
# with devtools::test() and on GitHub Actions.
skip_on_cran()

## helper: run the finalizers and retry pending deletions
collect <- function() {
  invisible(gc(verbose = FALSE)); invisible(gc(verbose = FALSE))
  dartR.base:::.fbm_retry_pending()
}

test_that("backing files of discarded subsets are deleted, the live one kept", {
  x <- gl.gen2fbm(testset.gl, verbose = 0)
  files <- x@fbm$backingfile
  for (k in 1:5) {
    x <- x[-1, ]
    files <- c(files, x@fbm$backingfile)
  }
  collect()
  expect_identical(file.exists(files), c(rep(FALSE, 5), TRUE))
  expect_identical(as.matrix(x), as.matrix(testset.gl[-(1:5), ]))
})

test_that("a file is kept while any copy of the object exists", {
  x <- gl.gen2fbm(testset.gl, verbose = 0)
  f <- x@fbm$backingfile
  y <- x
  y@other$note <- "modified copy"
  rm(x); collect()
  expect_true(file.exists(f))
  expect_identical(as.matrix(y), as.matrix(testset.gl))
  rm(y); collect()
  expect_false(file.exists(f))
})

test_that("a second R object on the same file keeps the file alive", {
  G <- utils.fbm.new(10, 10)
  f <- G$backingfile
  G2 <- dartR.base:::.fbm_register(
    bigsnpr::snp_fastImputeSimple(G, method = "mean2"))
  rm(G); collect()
  expect_true(file.exists(f))
  rm(G2); collect()
  expect_false(file.exists(f))
})

test_that("registering the same object twice does not keep its file forever", {
  G <- utils.fbm.new(5, 5)
  f <- G$backingfile
  dartR.base:::.fbm_register(G)
  dartR.base:::.fbm_register(G)
  rm(G); collect()
  expect_false(file.exists(f))
})

test_that("a backing file named by the user is never deleted", {
  stem <- file.path(tempdir(), "keep_me_genotypes")
  x <- gl.gen2fbm(testset.gl, backingfile = stem, verbose = 0)
  f <- x@fbm$backingfile
  rm(x); collect()
  expect_true(file.exists(f))
  unlink(f)
})

test_that("gl.fbm.cleanup deletes unused files and reports", {
  x <- gl.gen2fbm(testset.gl, verbose = 0)
  for (k in 1:3) x <- x[-1, ]
  res <- gl.fbm.cleanup(verbose = 0)
  expect_type(res, "list")
  expect_named(res, c("deleted", "deleted.mb", "in.use", "in.use.mb",
                      "pending"))
  expect_true(file.exists(x@fbm$backingfile))
  expect_equal(res$pending, 0)
  expect_gte(res$in.use, 1)
})

test_that("gl.pcoa on an FBM with missing data leaves no temporary file", {
  x <- gl.gen2fbm(testset.gl, verbose = 0)
  d <- getOption("FBM.dir", tempdir())
  before <- list.files(d, "^gl_pcoa_.*\\.bk$")
  p <- gl.pcoa(x, verbose = 0)
  collect()
  expect_identical(list.files(d, "^gl_pcoa_.*\\.bk$"), before)
})
