# Characterization tests for gl.rename.pop
# Baseline snapshotted before review (review-gl.rename.pop).
# These tests capture what the function DOES; assertions marked [approved diff]
# were flipped in Phase C to reflect approved behaviour changes.

test_that("gl.rename.pop renames a population exactly", {
  res <- gl.rename.pop(testset.gl, old = "EmsubRopeMata", new = "Outgroup",
                       verbose = 0)
  expect_true("Outgroup" %in% popNames(res))
  expect_false("EmsubRopeMata" %in% popNames(res))
  expect_equal(nPop(res), nPop(testset.gl))
  expect_equal(nInd(res), nInd(testset.gl))
  expect_identical(
    indNames(res)[as.character(pop(res)) == "Outgroup"],
    indNames(testset.gl)[as.character(pop(testset.gl)) == "EmsubRopeMata"]
  )
})

test_that("gl.rename.pop appends one history entry and leaves input untouched", {
  h <- length(testset.gl@other$history)
  res <- gl.rename.pop(testset.gl, old = "EmsubRopeMata", new = "Outgroup",
                       verbose = 0)
  expect_equal(length(res@other$history), h + 1)
  xcopy <- testset.gl
  invisible(gl.rename.pop(xcopy, old = "EmsubRopeMata", new = "Outgroup",
                          verbose = 0))
  expect_identical(xcopy, testset.gl)
})

test_that("gl.rename.pop works on SilicoDArT", {
  res <- gl.rename.pop(testset.gs, old = popNames(testset.gs)[1],
                       new = "Renamed1", verbose = 0)
  expect_true("Renamed1" %in% popNames(res))
})

test_that("gl.rename.pop missing old/new arguments error", {
  expect_error(gl.rename.pop(testset.gl, new = "X", verbose = 0))
  expect_error(gl.rename.pop(testset.gl, old = "EmmacRoss", verbose = 0))
})

test_that("gl.rename.pop with a nonexistent old population", {
  # [approved diff F1] a nonexistent population is now a fatal error
  # (previously a silent no-op that still appended a history entry).
  expect_error(
    gl.rename.pop(testset.gl, old = "NoSuchPop", new = "X", verbose = 0),
    "not present"
  )
})

test_that("gl.rename.pop renaming to an existing population", {
  # [approved diff F2, decision: forbid] renaming to an existing name is
  # now a fatal error (previously silently merged the two populations);
  # amalgamation goes through gl.recode.pop.
  expect_error(
    gl.rename.pop(testset.gl, old = "EmmacRoss", new = "EmmacTweeUki",
                  verbose = 0),
    "already exists"
  )
})

test_that("gl.rename.pop on a pop-less object", {
  # [approved diff F3] clear fatal error instead of 'attempt to set an
  # attribute on NULL'.
  xnp <- testset.gl
  xnp@pop <- NULL
  expect_error(
    gl.rename.pop(xnp, old = "A", new = "B", verbose = 0),
    "Population names not detected"
  )
})

test_that("gl.rename.pop returns invisibly", {
  # [approved diff F5] now returns invisibly.
  v <- withVisible(gl.rename.pop(testset.gl, old = "EmsubRopeMata",
                                 new = "Outgroup", verbose = 0))
  expect_false(v$visible)
})
