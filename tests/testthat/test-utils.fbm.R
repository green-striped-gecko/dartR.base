test_that("utils.fbm.tempfile honours option FBM.dir and creates the folder", {
  d <- file.path(tempdir(), "fbm_dir_test", "sub")
  unlink(dirname(d), recursive = TRUE)
  old <- options(FBM.dir = d)
  on.exit(options(old), add = TRUE)
  f <- utils.fbm.tempfile("abc_")
  expect_true(dir.exists(d))
  expect_identical(normalizePath(dirname(f)), normalizePath(d))
  expect_match(basename(f), "^abc_")
})

test_that("FBM objects made by dartR land in FBM.dir and are registered", {
  d <- file.path(tempdir(), "fbm_dir_test2")
  old <- options(FBM.dir = d)
  on.exit(options(old), add = TRUE)
  x <- gl.gen2fbm(testset.gl, verbose = 0)
  expect_identical(normalizePath(dirname(x@fbm$backingfile)), normalizePath(d))
  expect_true(dartR.base:::.fbm_is_registered(x@fbm))

  # subset, cbind, rbind, impute copy all go through the helpers
  y <- x[1:10, 1:20]
  expect_identical(normalizePath(dirname(y@fbm$backingfile)), normalizePath(d))
  expect_true(dartR.base:::.fbm_is_registered(y@fbm))
  cb <- cbind(x[, 1:5], x[, 6:10])
  expect_true(dartR.base:::.fbm_is_registered(cb@fbm))
  rb <- rbind(x[1:5, ], x[6:10, ])
  expect_identical(normalizePath(dirname(rb@fbm$backingfile)), normalizePath(d))
  expect_true(dartR.base:::.fbm_is_registered(rb@fbm))
})

test_that("a backing file named by the user is used and not registered", {
  stem <- file.path(tempdir(), "my_own_genotypes")
  x <- gl.gen2fbm(testset.gl, backingfile = stem, verbose = 0)
  expect_identical(normalizePath(x@fbm$backingfile),
                   normalizePath(paste0(stem, ".bk")))
  expect_false(dartR.base:::.fbm_is_registered(x@fbm))
})

test_that("rbind uses its backingfile argument", {
  x <- gl.gen2fbm(testset.gl, verbose = 0)
  stem <- file.path(tempdir(), "my_rbind")
  rb <- rbind(x[1:5, ], x[6:10, ], backingfile = stem)
  expect_identical(normalizePath(rb@fbm$backingfile),
                   normalizePath(paste0(stem, ".bk")))
})

test_that("utils.fbm.copy keeps genotypes and code", {
  x <- gl.gen2fbm(testset.gl, verbose = 0)
  G2 <- utils.fbm.copy(x@fbm, ind.row = 3:7, ind.col = 10:20)
  expect_identical(G2[], x@fbm[3:7, 10:20])
  expect_identical(G2$code256, x@fbm$code256)
  expect_false(identical(G2$backingfile, x@fbm$backingfile))
})
