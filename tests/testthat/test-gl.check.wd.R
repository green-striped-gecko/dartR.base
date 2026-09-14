# Characterization tests for gl.check.wd
# Baseline snapshotted before review (review-gl.check.wd).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.check.wd resolves the three documented sources", {
  td <- tempdir()
  old <- getOption("dartR_wd"); on.exit(options(dartR_wd = old), add = TRUE)
  options(dartR_wd = NULL)
  out <- capture.output(r <- gl.check.wd(verbose = 0))
  expect_equal(normalizePath(r), normalizePath(td))
  expect_length(out, 0)
  custom <- file.path(td, "cwd-src"); dir.create(custom, showWarnings = FALSE)
  options(dartR_wd = custom)
  r2 <- gl.check.wd(verbose = 0)
  expect_equal(normalizePath(r2), normalizePath(custom))
  options(dartR_wd = NULL)
  r3 <- gl.check.wd(custom, verbose = 0)
  expect_equal(normalizePath(r3), normalizePath(custom))
})

test_that("gl.check.wd falls back to tempdir for a non-existent path", {
  td <- tempdir()
  bad <- file.path(td, "no-such-dir-xyz")
  # [approved diff F1] the fallback warning is now gated at verbose >= 1
  out0 <- capture.output(r <- gl.check.wd(bad, verbose = 0))
  expect_equal(normalizePath(r), normalizePath(td))
  expect_length(out0, 0)
  out1 <- capture.output(r1 <- gl.check.wd(bad, verbose = 1))
  expect_true(any(grepl("does not exist", out1)))
})

test_that("gl.check.wd falls back to tempdir for a file (not a directory)", {
  td <- tempdir()
  f <- file.path(td, "cwd-file.txt"); writeLines("x", f)
  out <- capture.output(r <- gl.check.wd(f, verbose = 0))
  expect_equal(normalizePath(r), normalizePath(td))
})

test_that("gl.check.wd falls back to tempdir for a non-character wd", {
  # [approved diff F2] non-character wd no longer errors; it takes the
  # documented tempdir fallback
  td <- tempdir()
  out5 <- capture.output(r5 <- gl.check.wd(5, verbose = 0))
  expect_equal(normalizePath(r5), normalizePath(td))
  outna <- capture.output(rna <- gl.check.wd(NA, verbose = 0))
  expect_equal(normalizePath(rna), normalizePath(td))
})

test_that("gl.check.wd falls back to tempdir for a length>1 wd", {
  # [approved diff F3] multi-element wd no longer raises a condition-length
  # error; it takes the tempdir fallback
  td <- tempdir()
  custom <- file.path(td, "cwd-src2"); dir.create(custom, showWarnings = FALSE)
  out <- capture.output(r <- gl.check.wd(c(custom, td), verbose = 0))
  expect_equal(normalizePath(r), normalizePath(td))
})

test_that("gl.check.wd does not create a missing directory", {
  td <- tempdir()
  newdir <- file.path(td, "cwd-brand-new")
  if (dir.exists(newdir)) unlink(newdir, recursive = TRUE)
  out <- capture.output(r <- gl.check.wd(newdir, verbose = 0))
  expect_false(dir.exists(newdir))
  expect_equal(normalizePath(r), normalizePath(td))
})

test_that("gl.check.wd prints flag start/end and wd at verbose 2", {
  out <- capture.output(r <- gl.check.wd(verbose = 2))
  expect_true(any(grepl("Starting gl.check.wd", out)))
  expect_true(any(grepl("Working directory", out)))
  expect_true(any(grepl("Completed", out)))
})
