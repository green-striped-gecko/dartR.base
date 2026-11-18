# Characterization tests for utils.plot.save
# Baseline snapshotted before review (infrastructure wave, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.

test_that("saves an RDS and returns invisibly", {
  p <- ggplot2::ggplot()
  td <- file.path(tempdir(), "ps1")
  dir.create(td, showWarnings = FALSE)
  o <- capture.output(r <- utils.plot.save(p, dir = td, file = "t1",
                                           verbose = 0))
  expect_true(file.exists(file.path(td, "t1.RDS")))
  expect_null(r)
  expect_equal(length(o), 0)
})

test_that("the default verbose does not crash", {
  p <- ggplot2::ggplot()
  td <- file.path(tempdir(), "ps2")
  dir.create(td, showWarnings = FALSE)
  # [approved diff I2] baseline: verbose was never normalized, so the
  # default NULL crashed with "argument is of length zero".
  expect_error(utils.plot.save(p, dir = td, file = "t2"),
    NA)  # [approved diff I2]
})

test_that("a nonexistent directory falls back to a usable location", {
  p <- ggplot2::ggplot()
  # [approved diff I3] baseline: the fallback used tempfile() - a
  # nonexistent path - as the directory, so saveRDS crashed with
  # "cannot open the connection".
  expect_error(
    capture.output(utils.plot.save(p, dir = "Z:/no/such/dir9", file = "t3",
                                   verbose = 2)),
    NA)  # [approved diff I3]
})

test_that("file = NULL respects verbose 0", {
  p <- ggplot2::ggplot()
  # [approved diff I4] baseline: the "No plot saved" note was ungated.
  o <- capture.output(r <- utils.plot.save(p, file = NULL, verbose = 0))
  expect_equal(length(o), 0)  # [approved diff I4]
})
