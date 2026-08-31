test_that("gl.report.allna reports the all-NA loci of testset.gl correctly", {
  x <- testset.gl
  m <- as.matrix(x)
  expected <- locNames(x)[apply(m, 2, function(col) all(is.na(col)))]
  out <- capture.output(x2 <- gl.report.allna(x, verbose = 2))
  loc_line <- grep("loci that are missing", out, value = TRUE)
  expect_true(length(loc_line) >= 1)
  for (nm in expected) expect_true(any(grepl(nm, out, fixed = TRUE)))
})

test_that("gl.report.allna returns the input identical, no history append (SNP + SilicoDArT)", {
  x <- testset.gl
  out <- capture.output(x2 <- gl.report.allna(x, verbose = 0))
  expect_identical(x2, x)
  out_gs <- capture.output(x3 <- gl.report.allna(testset.gs, verbose = 0))
  expect_identical(x3, testset.gs)
})

test_that("gl.report.allna by.pop=TRUE per-population counts match independent computation", {
  x <- testset.gl
  out <- capture.output(x2 <- gl.report.allna(x, by.pop = TRUE, verbose = 2))
  # independent: union of loci all-NA within any one population
  expected <- unique(unlist(lapply(levels(pop(x)), function(p) {
    sub <- as.matrix(x[pop(x) == p, ])
    locNames(x)[apply(sub, 2, function(col) all(is.na(col)))]
  })))
  total_line <- grep("Loci all NA in one or more populations:", out, value = TRUE)
  expect_length(total_line, 1)
  expect_equal(as.numeric(sub(".*populations:\\s*", "", total_line)), length(expected))
})

test_that("verbose = 0 is fully silent, both branches (F3 fix)", {
  x <- testset.gl
  out0 <- capture.output(gl.report.allna(x, verbose = 0))
  expect_length(out0, 0)
  out0p <- capture.output(gl.report.allna(x, by.pop = TRUE, verbose = 0))
  expect_length(out0p, 0)
})

test_that("all-NA individuals listing names only the affected individuals, no NULLs (F1 fix)", {
  m2 <- matrix(c(0, 1, 2, 0, NA, NA, NA, NA, 2, 1, 0, 1), nrow = 3, byrow = TRUE,
               dimnames = list(c("IND_A", "IND_ALLNA", "IND_C"), paste0("L", 1:4)))
  gl_small <- new("genlight", gen = m2, ploidy = 2)
  pop(gl_small) <- factor(rep("p1", nInd(gl_small)))
  gl_small <- gl.compliance.check(gl_small, verbose = 0)
  outb <- capture.output(gl.report.allna(gl_small, verbose = 2))
  ind_line <- grep("individuals that are missing", outb, value = TRUE)
  expect_length(ind_line, 1)
  expect_true(grepl("IND_ALLNA", ind_line))
  expect_false(any(grepl("NULL", outb)))
})

test_that("non-genlight input fails fast with the standard message (F2 fix)", {
  expect_error(
    capture.output(gl.report.allna(data.frame(a = 1), verbose = 0)),
    "inappropriate object passed"
  )
})
