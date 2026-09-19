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

# Tier 3: depth pattern. A contaminant's alleles are a minority of the reads
# and are called mainly at loci with enough depth; a hybrid's alleles sit at
# 50 % and are called at any depth. No packaged dataset has both rdepth and
# two populations with enough fixed differences, so build two populations of
# 20 with 600 fixed differences and a lognormal depth, and give the first
# animal of population A the B allele at the fixed loci, either with a
# probability that rises with depth (contaminant) or at random (hybrid).
depth.host <- function(prob.fun, seed = 1) {
  set.seed(seed)
  n.a <- 20; n.b <- 20; n.loc <- 2000; n.diag <- 600
  p.a <- runif(n.loc, 0.05, 0.5); p.b <- runif(n.loc, 0.05, 0.5)
  g.a <- matrix(rbinom(n.a * n.loc, 2, rep(p.a, each = n.a)), n.a, n.loc)
  g.b <- matrix(rbinom(n.b * n.loc, 2, rep(p.b, each = n.b)), n.b, n.loc)
  g.a[, 1:n.diag] <- 0; g.b[, 1:n.diag] <- 2
  m <- rbind(g.a, g.b)
  depth <- exp(rnorm(n.loc, log(15), 0.6))
  hit <- which(runif(n.diag) < prob.fun(depth[1:n.diag]))
  m[1, hit] <- 1
  x <- new("genlight", m, ploidy = 2)
  indNames(x) <- paste0("i", seq_len(n.a + n.b))
  locNames(x) <- paste0("L", seq_len(n.loc))
  pop(x) <- rep(c("A", "B"), c(n.a, n.b))
  x@other$loc.metrics <- data.frame(rdepth = depth)
  list(x = x, host = "i1")
}

test_that("foreign alleles that rise with locus depth are reported as a dose pattern", {
  s <- depth.host(function(d) 0.6 * rank(d) / length(d))
  capture.output(r <- gl.report.contamination(s$x, plot.display = FALSE,
                                              verbose = 0))
  h <- r$ind[r$ind$id == s$host, ]
  expect_true(h$flag %in% c("suspect", "adjacent"))
  expect_equal(h$partner.pop, "B")
  expect_gt(h$depth.ratio, 1.5)
  expect_lt(h$depth.p, 0.01)
  expect_equal(h$pattern, "dose")
})

test_that("foreign alleles independent of locus depth are reported as a flat pattern", {
  s <- depth.host(function(d) rep(0.3, length(d)))
  capture.output(r <- gl.report.contamination(s$x, plot.display = FALSE,
                                              verbose = 0))
  h <- r$ind[r$ind$id == s$host, ]
  expect_true(h$flag %in% c("suspect", "adjacent"))
  expect_lt(h$depth.ratio, 1.25)
  expect_equal(h$pattern, "flat")
  # unflagged individuals carry no pattern
  expect_true(all(is.na(r$ind$pattern[r$ind$flag == ""])))
})

test_that("the depth pattern is NA when the genlight has no read depth", {
  x <- bandicoot.gl
  expect_null(x@other$loc.metrics$rdepth)
  out <- capture.output(r <- gl.report.contamination(x, plot.display = FALSE, verbose = 2))
  expect_true(all(is.na(r$ind$pattern)))
  expect_true(all(is.na(r$ind$depth.ratio)))
  expect_true(any(grepl("rdepth", out)))
})

test_that("plate adjacency is keyed by service when the column exists", {
  x <- platypus.gl
  im <- x@other$ind.metrics
  im$service <- NA; im$plate_location <- NA
  im$service[1:4] <- c("S1", "S1", "S2", "S2")
  im$plate_location[1:4] <- c("1-A1", "1-A2", "1-A1", "1-A2")
  x@other$ind.metrics <- im
  capture.output(r <- gl.report.contamination(x, plot.display = FALSE, verbose = 0))
  # A1 and A2 are adjacent within each order; "1-A1" of S1 is not adjacent
  # to "1-A2" of S2
  expect_equal(nrow(r$pairs), 2)
  expect_setequal(paste(r$pairs$id1, r$pairs$id2),
                  c(paste(indNames(x)[1], indNames(x)[2]),
                    paste(indNames(x)[3], indNames(x)[4])))
})

test_that("a second, heavily contaminated animal in the host population does not hide the dose pattern", {
  # Give i2 of population A the B allele at 90 % of the fixed loci (a heavy
  # contamination) before injecting the dose pattern into i1. Reference sets
  # built from unflagged members keep the diagnostic loci; leave-one-out
  # sets lose every locus at which i2 carries the B allele.
  s <- depth.host(function(d) 0.6 * rank(d) / length(d))
  m <- as.matrix(s$x)
  m[2, sample(1:600, 540)] <- 1
  x <- s$x
  x@gen <- new("genlight", m, ploidy = 2)@gen
  indNames(x) <- indNames(s$x); pop(x) <- pop(s$x)
  capture.output(r <- gl.report.contamination(x, plot.display = FALSE,
                                              verbose = 0))
  h <- r$ind[r$ind$id == s$host, ]
  expect_equal(h$pattern, "dose")
  expect_gt(h$foreign.rate, 0.2)
  expect_lt(h$foreign.rate, 0.5)
})

test_that("a mixture too heavy for depth to limit the calls is reported as saturated", {
  s <- depth.host(function(d) rep(0.9, length(d)))
  capture.output(r <- gl.report.contamination(s$x, plot.display = FALSE,
                                              verbose = 0))
  h <- r$ind[r$ind$id == s$host, ]
  expect_gt(h$foreign.rate, 0.8)
  expect_equal(h$pattern, "saturated")
})
