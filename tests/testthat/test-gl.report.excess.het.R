# Characterization baseline (function-review, 2026-09-24): outputs of the
# reviewed state. Detects change; does not assert correctness.
lbp.excess <- c("28680077-13-C/T", "28681088-13-A/G", "28685997-50-A/G",
                "28687640-57-C/T", "28691062-7-C/T", "28691064-10-A/G")

test_that("baseline: gl.report.excess.het loci on LBP and testset.gl", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  for (Y in c(TRUE, FALSE)) {
    capture.output(r <- suppressWarnings(
      gl.report.excess.het(LBP, Yates = Y, plot.display = FALSE, verbose = 0)))
    expect_setequal(r$removed.loci, lbp.excess)
    expect_equal(nrow(r$results.table), 7)
  }
  capture.output(r <- suppressWarnings(
    gl.report.excess.het(testset.gl, plot.display = FALSE, verbose = 0)))
  expect_setequal(r$removed.loci, c("100049990-20-G/T", "100050106-50-T/G",
                                    "100050129-31-G/A"))
  expect_equal(nrow(r$results.table), 33)
})

test_that("baseline: gl.report.excess.het leaves its input untouched", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  x <- LBP
  capture.output(suppressWarnings(
    gl.report.excess.het(x, plot.display = FALSE, verbose = 0)))
  expect_identical(x, LBP)
})

test_that("the deprecation advice reproduces the wrapper's loci", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  for (Y in c(TRUE, FALSE)) {
    msg <- NULL
    capture.output(r <- withCallingHandlers(
      gl.report.excess.het(LBP, Yates = Y, plot.display = FALSE, verbose = 0),
      warning = function(w) {
        if (grepl("deprecated", conditionMessage(w))) msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }))
    advised <- sub("^Use (.*) instead\\.$", "\\1", strsplit(msg, "\n")[[1]][2])
    advised <- sub("gl.report.hwe(", "gl.report.hwe(LBP, plot.out = FALSE, verbose = 0, ",
                   advised, fixed = TRUE)
    capture.output(r2 <- suppressWarnings(eval(parse(text = advised))))
    expect_setequal(unique(as.character(r2$Locus)), r$removed.loci)
  }
})

test_that("ignored plot arguments are named in a warning", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  out <- capture.output(suppressWarnings(
    gl.report.excess.het(LBP, plot.display = FALSE, plot.file = "p",
                         plot.dir = tempdir(), verbose = 1)))
  expect_true(any(grepl("ignored.*plot.file, plot.dir", out)))
  out <- capture.output(suppressWarnings(
    gl.report.excess.het(LBP, plot.display = FALSE, verbose = 1)))
  expect_false(any(grepl("ignored", out)))
})
