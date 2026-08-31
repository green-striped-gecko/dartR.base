test_that("gl.report.bases returns the input object unchanged (SNP)", {
  x <- testset.gl
  out <- capture.output(x2 <- gl.report.bases(x, plot.display = FALSE, verbose = 0))
  expect_identical(x2, x)
  expect_length(x2@other$history, length(x@other$history))
})

test_that("gl.report.bases returns the input object unchanged (SilicoDArT)", {
  x <- testset.gs
  out <- capture.output(x2 <- gl.report.bases(x, plot.display = FALSE, verbose = 0))
  expect_identical(x2, x)
})

test_that("gl.report.bases errors clearly when TrimmedSequence is absent", {
  x <- testset.gl
  x@other$loc.metrics$TrimmedSequence <- NULL
  expect_error(
    gl.report.bases(x, plot.display = FALSE, verbose = 0),
    "TrimmedSequence"
  )
})

test_that("gl.report.bases base frequencies match an independent computation (SNP)", {
  x <- testset.gl
  seqs <- as.character(x@other$loc.metrics$TrimmedSequence)
  total <- sum(nchar(seqs))
  A_pct <- sum(lengths(regmatches(seqs, gregexpr("A", seqs, fixed = TRUE)))) * 100 / total
  out <- capture.output(gl.report.bases(x, plot.display = FALSE, verbose = 1))
  a_line <- grep("^\\s+A:", out, value = TRUE)
  expect_length(a_line, 1)
  reported_A <- as.numeric(sub(".*A:\\s*", "", a_line))
  expect_equal(reported_A, round(A_pct, 2))
})

test_that("gl.report.bases ts/tv reported for SNP, not for SilicoDArT", {
  out_snp <- capture.output(gl.report.bases(testset.gl, plot.display = FALSE, verbose = 1))
  expect_true(any(grepl("Transitions", out_snp)))
  out_gs <- capture.output(gl.report.bases(testset.gs, plot.display = FALSE, verbose = 1))
  expect_false(any(grepl("Transitions", out_gs)))
})

test_that("verbose = 0 is fully silent; verbose >= 1 prints the results table (F1 fix, VRB5)", {
  out0 <- capture.output(gl.report.bases(testset.gl, plot.display = FALSE, verbose = 0))
  expect_length(out0, 0)
  out1 <- capture.output(gl.report.bases(testset.gl, plot.display = FALSE, verbose = 1))
  expect_true(any(grepl("Base frequencies", out1)))
})
