# Characterization tests for gl.join
# Baseline snapshotted before review (review-gl.join).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("join by shared loci (different individuals)", {
  x1 <- testset.gl[1:7, ]; x2 <- testset.gl[11:14, ]
  o <- capture.output(g <- gl.join(x1, x2, verbose = 0))
  expect_length(o, 0)
  expect_equal(nInd(g), 11)
  expect_equal(nLoc(g), nLoc(testset.gl))
  # [approved diff F1] baseline: the combined object LOST ind.metrics
  # (plain rbind returns NULL and gl.join never rebuilt it).
  expect_equal(nrow(g@other$ind.metrics), 11)
  expect_equal(as.character(g@other$ind.metrics$id), indNames(g))
})

test_that("join by shared individuals (different loci)", {
  y1 <- platypus.gl[, 1:100]; y2 <- platypus.gl[, 101:200]
  o <- capture.output(g <- gl.join(y1, y2, verbose = 0))
  expect_length(o, 0)
  expect_equal(nLoc(g), 200)
  expect_equal(nInd(g), nInd(platypus.gl))
  expect_equal(nrow(g@other$loc.metrics), 200)
})

test_that("mismatched objects are fatal", {
  expect_error(capture.output(gl.join(testset.gl[1:7, ],
                                      platypus.gl[1:4, ], verbose = 0)))
})

test_that("legacy method parameter", {
  # [approved diff F3] baseline: the deprecation warning printed at
  # verbose 0.
  x1 <- testset.gl[1:7, ]; x2 <- testset.gl[11:14, ]
  o <- capture.output(g <- gl.join(x1, x2, method = "join.by.loc",
                                   verbose = 0))
  expect_length(o, 0)
  expect_equal(nInd(g), 11)
  expect_error(capture.output(gl.join(x1, x2, method = "bogus",
                                      verbose = 0)))
})

test_that("historical end2end/sidebyside legacy values are honoured", {
  # [amendment, member-nominated re-review 2026-08-31] the pre-refactor
  # gl.join accepted method='end2end'/'sidebyside' and the docs kept
  # describing them, but the legacy shim only mapped join.by.loc/ind -
  # real callers (dartR.popgen gl.assign.*) crashed. Now mapped, and the
  # requested join is validated against the data.
  x1 <- testset.gl[1:7, ]; x2 <- testset.gl[11:14, ]
  invisible(capture.output(g1 <- gl.join(x1, x2, method = "end2end",
                                         verbose = 0)))
  invisible(capture.output(g0 <- gl.join(x1, x2, verbose = 0)))
  expect_equal(nInd(g1), 11)
  expect_identical(as.matrix(g1), as.matrix(g0))
  y1 <- platypus.gl[, 1:100]; y2 <- platypus.gl[, 101:200]
  invisible(capture.output(g2 <- gl.join(y1, y2, method = "sidebyside",
                                         verbose = 0)))
  expect_equal(nLoc(g2), 200)
  # a requested join that contradicts the data fails clearly
  expect_error(capture.output(gl.join(x1, x2, method = "sidebyside",
                                      verbose = 0)),
               "do not match")
  expect_error(capture.output(gl.join(y1, y2, method = "end2end",
                                      verbose = 0)),
               "do not match")
})

test_that("ind-join on testset-style flags (no OneRatio/PIC columns)", {
  # [approved diff F2] baseline: the flags block assigned OneRatio/PIC
  # products that do not exist in testset.gl's flags data.frame ->
  # "replacement has 0 rows" crash. Now only shared flags are combined.
  z1 <- testset.gl[, 1:100]; z2 <- testset.gl[, 101:200]
  invisible(capture.output(g <- gl.join(z1, z2, verbose = 0)))
  expect_equal(nLoc(g), 200)
  expect_equal(nInd(g), nInd(testset.gl))
})

test_that("SNP + SilicoDArT join is fatal", {
  # [approved diff F4] datatypes were checked individually but never
  # compared.
  expect_error(capture.output(gl.join(testset.gl[1:7, ], testset.gs[1:7, ],
                                      verbose = 0)),
               "different datatypes")
})

test_that("overlapping individual names are made unique in a loc-join", {
  w1 <- testset.gl[1:5, ]; w2 <- testset.gl[3:8, ]
  invisible(capture.output(g <- gl.join(w1, w2, verbose = 0)))
  expect_equal(nInd(g), 11)
  expect_equal(length(unique(indNames(g))), 11)
})

test_that("history is carried and appended", {
  # Docs claim the history is cleared; actual behaviour inherits x1's
  # history and appends the gl.join call (documented in the review).
  y1 <- platypus.gl[, 1:100]; y2 <- platypus.gl[, 101:200]
  invisible(capture.output(g <- gl.join(y1, y2, verbose = 0)))
  h <- g@other$history[[length(g@other$history)]]
  expect_true(is.call(h))
  expect_gte(length(g@other$history), 1)
})
