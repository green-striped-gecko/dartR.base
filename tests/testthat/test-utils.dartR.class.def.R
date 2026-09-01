# Characterization tests for utils.dartR.class.def (dartR S4 class layer)
# Baseline snapshotted before review (branch review-utils.dartR.class.def,
# dev at ddaed27). Assertions marked [approved diff] were flipped in
# Phase C.

make_sub <- function() gl.compliance.check(testset.gl[1:10, 1:30], verbose = 0)

test_that("subset, rbind and cbind round-trip exactly (gen-backed)", {
  sub <- make_sub()
  x2 <- sub[1:5, 1:10]
  expect_equal(c(nInd(x2), nLoc(x2)), c(5, 10))
  expect_equal(nrow(x2@other$loc.metrics), 10)
  r1 <- rbind(sub[1:5, ], sub[6:10, ])
  expect_identical(unname(as.matrix(r1)), unname(as.matrix(sub)))
  c1 <- cbind(sub[, 1:15], sub[, 16:30])
  expect_identical(unname(as.matrix(c1)), unname(as.matrix(sub)))
})

test_that("glSum and glMean match adegenet on gen-backed objects", {
  sub <- make_sub()
  gl_ref <- sub; class(gl_ref) <- "genlight"
  expect_equal(glSum(sub), adegenet::glSum(gl_ref))
  expect_equal(glMean(sub), adegenet::glMean(gl_ref))
})

test_that("show displays loc.metrics detail regardless of ind.metrics", {
  sub <- make_sub()
  s1 <- sub
  s1@other$ind.metrics <- NULL
  # [approved diff G1] baseline: the loc.metrics detail line was gated on
  # ind.metrics being present (copy-paste in the show method).
  o <- capture.output(show(s1))
  expect_true(any(grepl("^    @other[$]loc.metrics:", o)))  # [approved diff G1]
  o2 <- capture.output(show(sub))
  expect_true(any(grepl("^    @other[$]loc.metrics:", o2)))
})

test_that("character j with an unmatched locus name fails informatively", {
  sub <- make_sub()
  # [approved diff G2] baseline: nomatch = 0 indices reached the SNPbin
  # subsetter, crashing with "Cannot subset a SNPbin with mixed
  # subscripts."
  expect_error(sub[, c(locNames(sub)[1], "bogus_locus")],
    "bogus_locus")  # [approved diff G2]
  x <- sub[, locNames(sub)[2:4]]
  expect_equal(nLoc(x), 3)
})

test_that("the FBM layer matches the gen-backed reference", {
  sub <- make_sub()
  fb <- gl.gen2fbm(sub)
  expect_equal(unname(as.matrix(fb)), unname(as.matrix(sub)))
  expect_equal(unname(as.matrix(fb[1:5, 1:10])),
               unname(as.matrix(sub[1:5, 1:10])))
  expect_equal(unname(glSum(fb)), unname(glSum(sub)))
  expect_equal(unname(glMean(fb)), unname(glMean(sub)))
  expect_equal(sum(lengths(NA.posi(fb))),
               sum(lengths(lapply(sub@gen, function(e) e@NA.posi))))
})

test_that("negative indices work on FBM-backed objects", {
  sub <- make_sub()
  fb <- gl.gen2fbm(sub)
  # [approved diff G3] baseline: negative indices reached big_copy
  # unnormalized and crashed.
  fneg <- fb[-1, ]  # [approved diff G3]
  expect_equal(nInd(fneg), 9)
  expect_equal(nInd(sub[-1, ]), 9)  # gen-backed negative indexing works
})
