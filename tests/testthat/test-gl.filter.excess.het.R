# Characterization baseline (function-review, 2026-09-24): outputs of the
# reviewed state. Detects change; does not assert correctness.
test_that("baseline: gl.filter.excess.het removals on LBP and testset.gl", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  lbp.excess <- c("28680077-13-C/T", "28681088-13-A/G", "28685997-50-A/G",
                  "28687640-57-C/T", "28691062-7-C/T", "28691064-10-A/G")
  for (Y in c(TRUE, FALSE)) {
    capture.output(f <- suppressWarnings(
      gl.filter.excess.het(LBP, Yates = Y, verbose = 0)))
    expect_setequal(setdiff(locNames(LBP), locNames(f)), lbp.excess)
    expect_equal(nrow(f@other$loc.metrics), nLoc(f))
  }
  capture.output(f <- suppressWarnings(
    gl.filter.excess.het(testset.gl, verbose = 0)))
  expect_equal(nLoc(f), 252)
})

test_that("the history records the wrapper's own call, which replays", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  capture.output(f <- suppressWarnings(
    gl.filter.excess.het(LBP, Yates = TRUE, verbose = 0)))
  expect_length(f@other$history, length(LBP@other$history) + 1)
  h <- f@other$history[[length(f@other$history)]]
  expect_equal(as.character(h[[1]]), "gl.filter.excess.het")
  capture.output(f2 <- suppressWarnings(eval(h)))
  expect_equal(locNames(f2), locNames(f))
})

test_that("the deprecation advice reproduces the wrapper's removals", {
  skip_if_not_installed("HardyWeinberg")
  pdf(NULL); on.exit(dev.off())
  for (Y in c(TRUE, FALSE)) {
    msg <- NULL
    capture.output(f <- withCallingHandlers(
      gl.filter.excess.het(LBP, Yates = Y, verbose = 0),
      warning = function(w) {
        if (grepl("deprecated", conditionMessage(w))) msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }))
    advised <- sub("^Use (.*) instead\\.$", "\\1", strsplit(msg, "\n")[[1]][2])
    advised <- sub("gl.filter.hwe(", "gl.filter.hwe(LBP, verbose = 0, ", advised,
                   fixed = TRUE)
    capture.output(f2 <- suppressWarnings(eval(parse(text = advised))))
    expect_equal(locNames(f2), locNames(f))
  }
})
