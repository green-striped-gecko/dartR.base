# Characterization tests for gl.sort
# Baseline snapshotted before review (review-gl.sort).
# Assertions marked [approved diff] were flipped in Phase C.

test_that("default pop sort orders populations alphabetically", {
  o <- capture.output(g <- gl.sort(bandicoot.gl, verbose = 0))
  expect_length(o, 0)
  expect_false(is.unsorted(as.character(pop(g))))
  expect_equal(nInd(g), nInd(bandicoot.gl))
  expect_setequal(indNames(g), indNames(bandicoot.gl))
  # ind.metrics travel with the sort (same permutation as the genotypes;
  # note bandicoot.gl's id column does not match indNames even pre-sort)
  idx <- order(pop(bandicoot.gl))
  expect_equal(as.character(g@other$ind.metrics$id),
               as.character(bandicoot.gl@other$ind.metrics$id[idx]))
})

test_that("pop sort with explicit order.by honors the given order", {
  ord <- c("WA", "SA", "VIC", "NSW", "QLD")
  invisible(capture.output(g <- gl.sort(bandicoot.gl, sort.by = "pop",
                                        order.by = ord, verbose = 0)))
  expect_equal(levels(pop(g)), ord)
  expect_equal(unique(as.character(pop(g))), ord)
  expect_error(capture.output(gl.sort(bandicoot.gl, sort.by = "pop",
        order.by = c("WA", "SA"), verbose = 0)))
})

test_that("ind sort by name and by external vector", {
  invisible(capture.output(g <- gl.sort(bandicoot.gl, sort.by = "ind",
                                        verbose = 0)))
  expect_false(is.unsorted(indNames(g)))
  miss <- rowSums(is.na(as.matrix(bandicoot.gl)))
  invisible(capture.output(g2 <- gl.sort(bandicoot.gl, sort.by = "ind",
                                         order.by = miss, verbose = 0)))
  expect_false(is.unsorted(rowSums(is.na(as.matrix(g2)))))
  expect_error(capture.output(gl.sort(bandicoot.gl, sort.by = "bogus",
                                      verbose = 0)))
})

test_that("history receives a call", {
  # [approved diff F1] baseline: the entry was appended as
  # c(match.call()), coercing the call to a list.
  invisible(capture.output(g <- gl.sort(bandicoot.gl, verbose = 0)))
  h <- g@other$history[[length(g@other$history)]]
  expect_true(is.call(h))
})

test_that("script end flag", {
  # [approved diff F2] baseline: no FLAG SCRIPT END block existed -
  # Completed: never printed at any verbosity.
  o <- capture.output(g <- gl.sort(bandicoot.gl, verbose = 2))
  expect_true(any(grepl("Completed", o)))
})
