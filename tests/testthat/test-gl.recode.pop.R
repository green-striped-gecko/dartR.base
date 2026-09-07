# Characterization tests for gl.recode.pop
# Baseline snapshotted before review (review-gl.recode.pop).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

recode_file <- function() {
  system.file("extdata", "testset_pop_recode.csv", package = "dartR.data")
}

delete_file <- function(del_pops = c("EmmacRoss", "EmmacTweeUki")) {
  rt <- read.csv(recode_file(), header = FALSE, stringsAsFactors = FALSE)
  rt[rt[, 1] %in% del_pops, 2] <- "Delete"
  f <- file.path(tempdir(), "recode_del.csv")
  write.table(rt, f, sep = ",", row.names = FALSE, col.names = FALSE,
              quote = FALSE)
  f
}

test_that("gl.recode.pop applies the recode mapping exactly", {
  rt <- read.csv(recode_file(), header = FALSE, stringsAsFactors = FALSE)
  map <- setNames(rt[, 2], rt[, 1])
  res <- gl.recode.pop(testset.gl, pop.recode = recode_file(), verbose = 0)
  expect_identical(as.character(pop(res)),
                   unname(map[as.character(pop(testset.gl))]))
  expect_equal(nInd(res), nInd(testset.gl))
  expect_equal(nLoc(res), nLoc(testset.gl))
})

test_that("gl.recode.pop deletes populations flagged Delete", {
  del_pops <- c("EmmacRoss", "EmmacTweeUki")
  expected_deleted <- indNames(testset.gl)[
    as.character(pop(testset.gl)) %in% del_pops]
  res <- gl.recode.pop(testset.gl, pop.recode = delete_file(), verbose = 0)
  expect_equal(nInd(res), nInd(testset.gl) - length(expected_deleted))
  expect_false(any(c("Delete", "delete") %in% popNames(res)))
  expect_false(any(expected_deleted %in% indNames(res)))
})

test_that("gl.recode.pop input untouched and error paths", {
  xcopy <- testset.gl
  invisible(capture.output(
    gl.recode.pop(xcopy, pop.recode = recode_file(), verbose = 0)
  ))
  expect_identical(xcopy, testset.gl)
  # recode table with a missing population errors
  rt <- read.csv(recode_file(), header = FALSE, stringsAsFactors = FALSE)
  f <- file.path(tempdir(), "recode_short.csv")
  write.table(rt[-1, ], f, sep = ",", row.names = FALSE, col.names = FALSE,
              quote = FALSE)
  expect_error(gl.recode.pop(testset.gl, pop.recode = f, verbose = 0))
})

test_that("gl.recode.pop history entries", {
  # [approved diff F3] one call now appends exactly one entry on both the
  # plain and Delete paths (pre-fix the Delete path appended two).
  h <- length(testset.gl@other$history)
  res1 <- gl.recode.pop(testset.gl, pop.recode = recode_file(), verbose = 0)
  expect_equal(length(res1@other$history), h + 1)
  res2 <- gl.recode.pop(testset.gl, pop.recode = delete_file(), verbose = 0)
  expect_equal(length(res2@other$history), h + 1)
  expect_match(deparse(res2@other$history[[h + 1]])[1], "gl.recode.pop")
})

test_that("gl.recode.pop deletions listing at verbose = 3", {
  # [approved diff F1] the listing now names the actual deletions with the
  # right count, and appears at verbose >= 3 (including 5).
  del_pops <- c("EmmacRoss", "EmmacTweeUki")
  expected_deleted <- indNames(testset.gl)[
    as.character(pop(testset.gl)) %in% del_pops]
  out <- capture.output(
    res <- gl.recode.pop(testset.gl, pop.recode = delete_file(), verbose = 3)
  )
  expect_true(all(sapply(expected_deleted[1:3], function(n)
    any(grepl(n, out, fixed = TRUE)))))
  expect_true(any(grepl("total of 20", out)))
  expect_false(any(grepl("total of 16", out)))
  out5 <- capture.output(
    res <- gl.recode.pop(testset.gl, pop.recode = delete_file(), verbose = 5)
  )
  expect_true(any(grepl("Dropping", out5)))
})

test_that("gl.recode.pop flags no longer depend on verbosity (plain recode)", {
  # [approved diff F2] a pure renaming run keeps its flags valid at every
  # verbosity (pre-fix, verbose >= 2 invalidated them).
  r0 <- gl.recode.pop(testset.gl, pop.recode = recode_file(), verbose = 0)
  expect_true(r0@other$loc.metrics.flags$CallRate)
  invisible(capture.output(
    r2 <- gl.recode.pop(testset.gl, pop.recode = recode_file(), verbose = 2)
  ))
  expect_true(r2@other$loc.metrics.flags$CallRate)
  # deletion runs still get their flags reset (by the gl.drop.pop
  # delegation), at every verbosity
  rdel <- gl.recode.pop(testset.gl, pop.recode = delete_file(), verbose = 0)
  expect_false(rdel@other$loc.metrics.flags$CallRate)
})

test_that("gl.recode.pop returns invisibly", {
  # [approved diff F4] now returns invisibly.
  v <- withVisible(gl.recode.pop(testset.gl, pop.recode = recode_file(),
                                 verbose = 0))
  expect_false(v$visible)
})
