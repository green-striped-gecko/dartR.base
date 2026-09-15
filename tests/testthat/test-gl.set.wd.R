# Characterization tests for gl.set.wd
# Baseline snapshotted before review (review-gl.set.wd).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect the approved behaviour change (error
# loudly on an invalid path).

reset_wd <- function() options(dartR_wd = NULL)

test_that("gl.set.wd sets the option and returns the path for a valid dir", {
  reset_wd(); on.exit(reset_wd(), add = TRUE)
  td <- tempdir()
  custom <- file.path(td, "setwd-ok"); dir.create(custom, showWarnings = FALSE)
  out <- capture.output(r <- gl.set.wd(custom, verbose = 0))
  expect_equal(normalizePath(r), normalizePath(custom))
  expect_equal(normalizePath(getOption("dartR_wd")), normalizePath(custom))
  expect_length(out, 0)
})

test_that("gl.set.wd default sets the option to a fresh tempdir", {
  reset_wd(); on.exit(reset_wd(), add = TRUE)
  td <- tempdir()
  r <- gl.set.wd(verbose = 0)
  expect_equal(normalizePath(r), normalizePath(td))
  expect_equal(normalizePath(getOption("dartR_wd")), normalizePath(td))
})

test_that("gl.set.wd errors on a non-existent path and leaves wd unchanged", {
  # [approved diff F1] previously returned the bad path and printed a false
  # success message while silently not setting the option
  reset_wd(); on.exit(reset_wd(), add = TRUE)
  td <- tempdir()
  prior <- file.path(td, "setwd-prior"); dir.create(prior, showWarnings = FALSE)
  options(dartR_wd = prior)
  bad <- file.path(td, "setwd-nonexistent")
  expect_error(gl.set.wd(bad, verbose = 0), "does not exist")
  # the prior wd is untouched
  expect_equal(normalizePath(getOption("dartR_wd")), normalizePath(prior))
})

test_that("gl.set.wd errors on a file (not a directory)", {
  # [approved diff F1]
  reset_wd(); on.exit(reset_wd(), add = TRUE)
  td <- tempdir()
  f <- file.path(td, "setwd-file.txt"); writeLines("x", f)
  expect_error(gl.set.wd(f, verbose = 0), "does not exist")
})

test_that("gl.set.wd errors clearly on wd = NULL and non-character wd", {
  # [approved diff F2] no longer an opaque "invalid filename argument"
  reset_wd(); on.exit(reset_wd(), add = TRUE)
  expect_error(gl.set.wd(wd = NULL, verbose = 0), "does not exist")
  expect_error(gl.set.wd(5, verbose = 0), "does not exist")
})

test_that("gl.set.wd errors clearly on a length>1 wd", {
  # [approved diff F3] no longer a condition-length error
  reset_wd(); on.exit(reset_wd(), add = TRUE)
  td <- tempdir()
  custom <- file.path(td, "setwd-v"); dir.create(custom, showWarnings = FALSE)
  expect_error(gl.set.wd(c(custom, td), verbose = 0), "does not exist")
})

test_that("gl.set.wd and gl.check.wd round-trip for a valid path", {
  reset_wd(); on.exit(reset_wd(), add = TRUE)
  td <- tempdir()
  custom <- file.path(td, "setwd-rt"); dir.create(custom, showWarnings = FALSE)
  gl.set.wd(custom, verbose = 0)
  r <- gl.check.wd(verbose = 0)
  expect_equal(normalizePath(r), normalizePath(custom))
})
