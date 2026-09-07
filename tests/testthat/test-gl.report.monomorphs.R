test_that("gl.report.monomorphs counts match an independent computation (SNP)", {
  x <- testset.gl
  m <- as.matrix(x)
  is_mono <- apply(m, 2, function(col) {
    v <- col[!is.na(col)]
    length(v) == 0 || all(v == 0) || all(v == 2)
  })
  all_na <- apply(m, 2, function(col) all(is.na(col)))
  out <- capture.output(x2 <- gl.report.monomorphs(x, verbose = 1))
  mono_line <- grep("Monomorphic loci:", out, value = TRUE)
  na_line <- grep("all NA:", out, value = TRUE)
  expect_equal(as.numeric(sub(".*Monomorphic loci:\\s*", "", mono_line)), sum(is_mono))
  expect_equal(as.numeric(sub(".*all NA:\\s*", "", na_line)), sum(all_na))
})

test_that("gl.report.monomorphs counts match an independent computation (SilicoDArT)", {
  x <- testset.gs
  m <- as.matrix(x)
  is_mono <- apply(m, 2, function(col) {
    v <- col[!is.na(col)]
    length(v) == 0 || all(v == 0) || all(v == 1)
  })
  out <- capture.output(x2 <- gl.report.monomorphs(x, verbose = 1))
  mono_line <- grep("Monomorphic loci:", out, value = TRUE)
  expect_equal(as.numeric(sub(".*Monomorphic loci:\\s*", "", mono_line)), sum(is_mono))
})

test_that("returned object is identical to the input; no history append (F1 fix)", {
  x <- testset.gl
  out <- capture.output(x2 <- gl.report.monomorphs(x, verbose = 1))
  expect_identical(x2, x)
  out_gs <- capture.output(x3 <- gl.report.monomorphs(testset.gs, verbose = 0))
  expect_identical(x3, testset.gs)
})

test_that("verbose = 0 is fully silent; verbose >= 1 prints the report (F2 fix)", {
  x <- testset.gl
  out0 <- capture.output(gl.report.monomorphs(x, verbose = 0))
  expect_length(out0, 0)
  out1 <- capture.output(gl.report.monomorphs(x, verbose = 1))
  expect_true(any(grepl("Monomorphic loci:", out1)))
})
