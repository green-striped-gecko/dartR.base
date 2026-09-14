# Characterization tests for gl.edit.recode.pop
# Baseline snapshotted before review (review-gl.edit.recode.pop).
# The interactive edit() step is mocked (in the dartR.base namespace, which
# imports utils::edit) to return its input unchanged, so the surrounding logic
# is exercised deterministically without opening an editor. Assertions marked
# [approved diff] were flipped in Phase C to reflect approved behaviour changes.

make_sub <- function(n = 30) {
  s <- testset.gl[1:n, ]
  s@other$loc.metrics <- testset.gl@other$loc.metrics
  s
}
edit_noop <- function(name = NULL, ...) name

test_that("gl.edit.recode.pop resets metric flags regardless of verbosity (recalc = FALSE)", {
  local_mocked_bindings(edit = edit_noop, .package = "dartR.base")
  # [approved diff F2] utils.reset.flags now runs whenever recalc = FALSE
  sub <- make_sub()
  expect_true(sub@other$loc.metrics.flags$CallRate)
  r0 <- gl.edit.recode.pop(sub, recalc = FALSE, verbose = 0)
  out <- capture.output(r2 <- gl.edit.recode.pop(sub, recalc = FALSE,
                                                 verbose = 2))
  expect_false(r0@other$loc.metrics.flags$CallRate)
  expect_false(r2@other$loc.metrics.flags$CallRate)
})

test_that("gl.edit.recode.pop tolerates a missing monomorphs flag", {
  local_mocked_bindings(edit = edit_noop, .package = "dartR.base")
  # [approved diff F3] no longer crashes on a NULL flag
  sub <- make_sub()
  sub@other$loc.metrics.flags$monomorphs <- NULL
  expect_no_error(capture.output(gl.edit.recode.pop(sub, verbose = 0)))
})

test_that("gl.edit.recode.pop default arguments are recalc = FALSE, mono.rm = FALSE", {
  fm <- formals(gl.edit.recode.pop)
  expect_false(eval(fm$recalc))
  expect_false(eval(fm$mono.rm))
})

test_that("gl.edit.recode.pop reads pop.recode as input and does not overwrite it", {
  local_mocked_bindings(edit = edit_noop, .package = "dartR.base")
  # [approved diff F1] pop.recode is now a read-only input, loaded and applied
  # (edit no-ops here, so the loaded table is applied verbatim); the input
  # file is not overwritten
  td <- tempdir()
  sub <- make_sub()
  levs <- levels(pop(sub))
  existing <- file.path(td, "existing_recode.csv")
  writeLines(paste0(levs, ",MERGED"), existing)
  before <- readLines(existing)
  out <- capture.output(
    r <- gl.edit.recode.pop(sub, pop.recode = existing, verbose = 0)
  )
  expect_true(all(as.character(pop(r)) == "MERGED"))
  expect_identical(readLines(existing), before)
})

test_that("gl.edit.recode.pop writes the edited table to out.recode.file under outpath", {
  local_mocked_bindings(edit = edit_noop, .package = "dartR.base")
  # [approved diff F1] out.recode.file is the output, written under outpath
  td <- tempdir()
  sub <- make_sub()
  subdir <- file.path(td, "prdir"); dir.create(subdir, showWarnings = FALSE)
  if (file.exists(file.path(subdir, "out.csv"))) file.remove(file.path(subdir, "out.csv"))
  owd <- getwd(); on.exit(setwd(owd), add = TRUE); setwd(td)
  out <- capture.output(
    gl.edit.recode.pop(sub, out.recode.file = "out.csv", outpath = subdir,
                       verbose = 0)
  )
  expect_true(file.exists(file.path(subdir, "out.csv")))
  expect_false(file.exists(file.path(td, "out.csv")))
})

test_that("gl.edit.recode.pop errors when x has no populations", {
  local_mocked_bindings(edit = edit_noop, .package = "dartR.base")
  sub <- make_sub()
  pop(sub) <- NULL
  expect_error(gl.edit.recode.pop(sub, verbose = 0),
               "Population names not detected")
})

test_that("gl.edit.recode.pop appends a history entry and returns a genlight", {
  local_mocked_bindings(edit = edit_noop, .package = "dartR.base")
  sub <- make_sub()
  h0 <- length(sub@other$history)
  out <- capture.output(r <- gl.edit.recode.pop(sub, verbose = 0))
  expect_s4_class(r, "genlight")
  expect_equal(length(r@other$history), h0 + 1)
})
