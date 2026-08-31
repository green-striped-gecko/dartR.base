# Characterization tests for gl.recode.ind
# Baseline snapshotted before review (review-gl.recode.ind).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

recode_file <- function() {
  system.file("extdata", "testset_ind_recode.csv", package = "dartR.data")
}

rename_file <- function() {
  rt <- read.csv(recode_file(), header = FALSE, stringsAsFactors = FALSE)
  rt[tolower(rt[, 2]) == "delete", 2] <- rt[tolower(rt[, 2]) == "delete", 1]
  rt[1, 2] <- "RENAMED_1"
  f <- file.path(tempdir(), "ind_recode_rename.csv")
  write.table(rt, f, sep = ",", row.names = FALSE, col.names = FALSE,
              quote = FALSE)
  f
}

test_that("gl.recode.ind applies the recode mapping exactly", {
  rt <- read.csv(recode_file(), header = FALSE, stringsAsFactors = FALSE)
  map <- setNames(rt[, 2], rt[, 1])
  out <- capture.output(
    res <- gl.recode.ind(testset.gl, ind.recode = recode_file(), verbose = 0)
  )
  expected <- unname(map[indNames(testset.gl)])
  expected <- expected[tolower(expected) != "delete"]
  expect_identical(indNames(res), expected)
  expect_equal(nInd(res), nInd(testset.gl) - 1)   # fixture deletes UC_00126c
  expect_false("UC_00126c" %in% indNames(res))
})

test_that("gl.recode.ind input untouched and short-table error", {
  xcopy <- testset.gl
  invisible(capture.output(
    gl.recode.ind(xcopy, ind.recode = recode_file(), verbose = 0)
  ))
  expect_identical(xcopy, testset.gl)
  rt <- read.csv(recode_file(), header = FALSE, stringsAsFactors = FALSE)
  f <- file.path(tempdir(), "ind_recode_short.csv")
  write.table(rt[-1, ], f, sep = ",", row.names = FALSE, col.names = FALSE,
              quote = FALSE)
  expect_error(gl.recode.ind(testset.gl, ind.recode = f, verbose = 0))
})

test_that("gl.recode.ind verbose = 0 output", {
  # [approved diff F1] verbose = 0 is now fully silent — only present
  # spellings are passed to gl.drop.ind, so its 'not present' warning no
  # longer fires.
  out <- capture.output(
    res <- gl.recode.ind(testset.gl, ind.recode = recode_file(), verbose = 0)
  )
  expect_length(out, 0)
})

test_that("gl.recode.ind history entries", {
  # [approved diff F3] one call now appends exactly one entry.
  h <- length(testset.gl@other$history)
  res <- gl.recode.ind(testset.gl, ind.recode = recode_file(), verbose = 0)
  expect_equal(length(res@other$history), h + 1)
  expect_match(deparse(res@other$history[[h + 1]])[1], "gl.recode.ind")
})

test_that("gl.recode.ind deletions listing at verbose = 3", {
  # [approved diff F1] the listing now names the deleted individual's
  # original identifier and appears at verbose >= 3 (including 5).
  out <- capture.output(
    res <- gl.recode.ind(testset.gl, ind.recode = recode_file(), verbose = 3)
  )
  drop_idx <- grep("Dropping", out)
  expect_true(any(grepl("UC_00126c", out[drop_idx + 1])))
  expect_true(any(grepl("total of 1", out)))
  out5 <- capture.output(
    res <- gl.recode.ind(testset.gl, ind.recode = recode_file(), verbose = 5)
  )
  expect_true(any(grepl("Dropping", out5)))
})

test_that("gl.recode.ind flags no longer depend on verbosity (pure rename)", {
  # [approved diff F2] flags stay valid at every verbosity on a pure
  # rename; deletion runs still get gl.drop.ind's reset.
  r0 <- gl.recode.ind(testset.gl, ind.recode = rename_file(), verbose = 0)
  expect_true(r0@other$loc.metrics.flags$CallRate)
  invisible(capture.output(
    r2 <- gl.recode.ind(testset.gl, ind.recode = rename_file(), verbose = 2)
  ))
  expect_true(r2@other$loc.metrics.flags$CallRate)
  rdel <- gl.recode.ind(testset.gl, ind.recode = recode_file(), verbose = 0)
  expect_false(rdel@other$loc.metrics.flags$CallRate)
})

test_that("gl.recode.ind summary gate", {
  # [approved diff F4] the results summary now gates at verbose >= 3,
  # aligned with gl.recode.pop.
  out <- capture.output(
    r <- gl.recode.ind(testset.gl, ind.recode = rename_file(), verbose = 2)
  )
  expect_false(any(grepl("Summary of recoded", out)))
  out3 <- capture.output(
    r <- gl.recode.ind(testset.gl, ind.recode = rename_file(), verbose = 3)
  )
  expect_true(any(grepl("Summary of recoded", out3)))
})

test_that("gl.recode.ind returns invisibly", {
  # [approved diff F5] now returns invisibly.
  v <- withVisible(gl.recode.ind(testset.gl, ind.recode = recode_file(),
                                 verbose = 0))
  expect_false(v$visible)
})
