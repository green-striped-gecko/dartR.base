# Characterization tests for gl.edit.recode.ind
# Baseline snapshotted before review (review-gl.edit.recode.ind).
# The interactive edit() step no-ops under a non-interactive session, so these
# tests exercise the surrounding logic (which runs regardless of edits).
# Assertions marked [approved diff] were flipped in Phase C to reflect approved
# behaviour changes.

make_sub <- function(n = 20) {
  s <- testset.gl[1:n, ]
  s@other$loc.metrics <- testset.gl@other$loc.metrics
  s
}

test_that("gl.edit.recode.ind writes the recode file to outpath", {
  # [approved diff F4] the file now goes to outfilespec (outpath), not getwd()
  td <- tempdir()
  sub <- make_sub()
  subdir <- file.path(td, "recodedir"); dir.create(subdir, showWarnings = FALSE)
  owd <- getwd(); on.exit(setwd(owd), add = TRUE)
  setwd(td)
  if (file.exists(file.path(td, "myrecode.csv"))) file.remove(file.path(td, "myrecode.csv"))
  if (file.exists(file.path(subdir, "myrecode.csv"))) file.remove(file.path(subdir, "myrecode.csv"))
  out <- capture.output(
    gl.edit.recode.ind(sub, out.recode.file = "myrecode.csv",
                       outpath = subdir, verbose = 0)
  )
  expect_true(file.exists(file.path(subdir, "myrecode.csv")))
  expect_false(file.exists(file.path(td, "myrecode.csv")))
})

test_that("gl.edit.recode.ind resets metric flags regardless of verbosity (recalc = FALSE)", {
  # [approved diff F1] utils.reset.flags now runs whenever recalc = FALSE,
  # so the flag state no longer depends on verbose
  sub <- make_sub()
  expect_true(sub@other$loc.metrics.flags$CallRate)
  r0 <- gl.edit.recode.ind(sub, recalc = FALSE, verbose = 0)
  out <- capture.output(r2 <- gl.edit.recode.ind(sub, recalc = FALSE,
                                                 verbose = 2))
  expect_false(r0@other$loc.metrics.flags$CallRate)
  expect_false(r2@other$loc.metrics.flags$CallRate)
})

test_that("gl.edit.recode.ind tolerates a missing monomorphs flag", {
  # [approved diff F2] no longer crashes on a NULL flag
  sub <- make_sub()
  sub@other$loc.metrics.flags$monomorphs <- NULL
  expect_no_error(
    capture.output(r <- gl.edit.recode.ind(sub, verbose = 0))
  )
})

test_that("gl.edit.recode.ind default arguments are recalc = FALSE, mono.rm = FALSE", {
  fm <- formals(gl.edit.recode.ind)
  expect_false(eval(fm$recalc))
  expect_false(eval(fm$mono.rm))
})

test_that("gl.edit.recode.ind appends a history entry and returns a genlight", {
  sub <- make_sub()
  h0 <- length(sub@other$history)
  out <- capture.output(r <- gl.edit.recode.ind(sub, verbose = 0))
  expect_s4_class(r, "genlight")
  expect_equal(length(r@other$history), h0 + 1)
})
