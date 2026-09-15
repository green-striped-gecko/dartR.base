test_that("gl.report.contamination returns the expected structure and is silent at verbose = 0", {
  x <- testset.gl
  out <- capture.output(r <- gl.report.contamination(x, plot.display = FALSE, verbose = 0))
  expect_length(out, 0)
  expect_type(r, "list")
  expect_named(r, c("ind", "pairs", "kinship"))
  expect_equal(nrow(r$ind), nInd(x))
  expect_setequal(r$ind$id, indNames(x))
  expect_true(all(r$ind$flag %in% c("", "suspect", "adjacent", "rare-only")))
  expect_equal(dim(r$kinship), c(nInd(x), nInd(x)))
  expect_true(all(is.na(diag(r$kinship))))
  # testset.gl carries no plate positions, so adjacency is not tested
  expect_null(r$pairs)
  expect_true(all(is.na(r$ind$adjacent)))
})

test_that("heterozygosity matches an independent computation", {
  x <- testset.gl
  capture.output(r <- gl.report.contamination(x, plot.display = FALSE, verbose = 0))
  m <- as.matrix(x)
  het <- rowMeans(m == 1, na.rm = TRUE)
  expect_equal(r$ind$het[match(indNames(x), r$ind$id)], unname(round(het, 4)))
})

test_that("a synthetic mixture of two individuals is flagged and its donor named", {
  # testset.gl has too few loci for this (host and donor differ at 4 loci);
  # bandicoot.gl has 1000 loci and five populations. The host is the WA animal
  # with median heterozygosity, the donor the first SA animal.
  set.seed(1)
  x <- bandicoot.gl
  m <- as.matrix(x)
  het0 <- rowMeans(m == 1, na.rm = TRUE)
  wa <- which(pop(x) == "WA")
  host <- wa[which.min(abs(het0[wa] - median(het0[wa])))]
  donor <- which(pop(x) == "SA")[1]
  capture.output(base <- gl.report.contamination(x, plot.display = FALSE, verbose = 0))
  base.sus <- base$ind$id[base$ind$flag %in% c("suspect", "adjacent")]
  expect_false(indNames(x)[host] %in% base.sus)
  # a 50% contaminant reads as heterozygous wherever host and donor differ
  loci <- which(!is.na(m[host, ]) & !is.na(m[donor, ]) & m[host, ] != m[donor, ])
  expect_gt(length(loci), 100)
  m[host, sample(loci, round(0.5 * length(loci)))] <- 1
  x2 <- x
  x2@gen <- new("genlight", m, ploidy = 2)@gen
  indNames(x2) <- indNames(x); pop(x2) <- pop(x)
  capture.output(r <- gl.report.contamination(x2, plot.display = FALSE, verbose = 0))
  h <- r$ind[r$ind$id == indNames(x)[host], ]
  expect_true(h$flag %in% c("suspect", "adjacent"))
  expect_gt(h$het.excess, 0.02)
  # tier 2 names the donor itself as the source
  expect_equal(h$top.partner, indNames(x)[donor])
  expect_equal(h$partner.pop, "SA")
  # the mixture adds no suspect other than the host (a borderline baseline
  # suspect may drop out because the host shifts its population's median)
  sus <- r$ind$id[r$ind$flag %in% c("suspect", "adjacent")]
  expect_equal(setdiff(sus, base.sus), indNames(x)[host])
})

test_that("plate adjacency is read from a supplied plate table", {
  x <- testset.gl
  ids <- indNames(x)
  wells <- paste0(rep(LETTERS[1:8], length.out = length(ids)),
                  rep(1:12, each = 8, length.out = length(ids)))
  plate <- data.frame(id = ids, plate = "1", well = wells)
  capture.output(r <- gl.report.contamination(x, plate = plate, plot.display = FALSE, verbose = 0))
  expect_false(is.null(r$pairs))
  expect_true(all(r$ind$adjacent %in% c(TRUE, FALSE)))
  expect_true(all(nchar(r$pairs$well1) >= 2))
})

test_that("bad parameters warn and coerce, or stop", {
  x <- testset.gl
  out <- capture.output(r <- gl.report.contamination(x, rare.freq = 0.9, plot.display = FALSE, verbose = 2))
  expect_true(any(grepl("rare.freq must lie", out)))
  expect_error(capture.output(gl.report.contamination(x, z.flag = 0, plot.display = FALSE, verbose = 0)))
  expect_error(capture.output(gl.report.contamination(x, min.excess = -1, plot.display = FALSE, verbose = 0)))
  expect_error(capture.output(gl.report.contamination(x, plate = data.frame(id = 1), plot.display = FALSE, verbose = 0)))
  expect_error(capture.output(gl.report.contamination(testset.gs, plot.display = FALSE, verbose = 0)))
})
