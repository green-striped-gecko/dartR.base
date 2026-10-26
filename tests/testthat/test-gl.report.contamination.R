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

test_that("bad parameters stop", {
  x <- testset.gl
  expect_error(capture.output(gl.report.contamination(x, rare.freq = 0.9, plot.display = FALSE, verbose = 0)),
               "rare.freq must lie")
  expect_error(capture.output(gl.report.contamination(x, min.n = 1, plot.display = FALSE, verbose = 0)),
               "min.n must be at least 2")
  expect_error(capture.output(gl.report.contamination(x, rare.freq = NA, plot.display = FALSE, verbose = 0)),
               "rare.freq must be a single finite number")
  expect_error(capture.output(gl.report.contamination(x, z.flag = c(3, 4), plot.display = FALSE, verbose = 0)),
               "z.flag must be a single finite number")
  expect_error(capture.output(gl.report.contamination(x, min.share = "20", plot.display = FALSE, verbose = 0)),
               "min.share must be a single finite number")
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

# Characterization baseline (function-review, 2026-09-24): outputs of the
# reviewed state on packaged data. Detects change; does not assert
# correctness.
test_that("baseline: flags and summary statistics on packaged data", {
  base <- list(
    testset.gl = list(id = c("AA019157", "AA011737"),
                      flag = c("suspect", "rare-only"),
                      het = 4.1433, rare = 0.4616, kin = 126.7559, nf = 97),
    platypus.gl = list(id = "SUS21", flag = "rare-only",
                       het = 10.1192, rare = 0.2994, kin = 11.7857, nf = 277),
    bandicoot.gl = list(id = c("bc82", "bc16", "bc38", "bc39", "bc43",
                               "bc74", "bc21"),
                        flag = rep("suspect", 7),
                        het = 30.9664, rare = 0.2872, kin = 3.2913, nf = 284))
  for (nm in names(base)) {
    b <- base[[nm]]
    capture.output(r <- gl.report.contamination(get(nm), plot.display = FALSE,
                                                verbose = 0))
    fl <- r$ind$flag != ""
    expect_equal(r$ind$id[fl], b$id, label = nm)
    expect_equal(r$ind$flag[fl], b$flag, label = nm)
    expect_equal(sum(r$ind$het, na.rm = TRUE), b$het, tolerance = 1e-4)
    expect_equal(sum(r$ind$rare.burden, na.rm = TRUE), b$rare,
                 tolerance = 1e-4)
    expect_equal(sum(r$kinship, na.rm = TRUE), b$kin, tolerance = 1e-4)
    expect_equal(sum(r$ind$n.foreign), b$nf)
  }
})

test_that("adjacent pairs match an all-pairs search", {
  x <- testset.gl
  ids <- indNames(x)
  n <- length(ids)
  # two plates, and two individuals sharing a well, to exercise the lookup
  wells <- paste0(LETTERS[(0:(n - 1)) %% 8 + 1], ((0:(n - 1)) %/% 8) %% 12 + 1)
  plate <- data.frame(id = ids, plate = (0:(n - 1)) %/% 96, well = wells)
  plate$well[2] <- plate$well[1]
  capture.output(r <- gl.report.contamination(x, plate = plate, plot.display = FALSE, verbose = 0))
  rw <- match(substr(plate$well, 1, 1), LETTERS)
  cl <- as.integer(substring(plate$well, 2))
  pr <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  adj <- plate$plate[pr[, 1]] == plate$plate[pr[, 2]] &
    abs(rw[pr[, 1]] - rw[pr[, 2]]) + abs(cl[pr[, 1]] - cl[pr[, 2]]) == 1
  expect_setequal(paste(r$pairs$id1, r$pairs$id2),
                  paste(ids[pr[adj, 1]], ids[pr[adj, 2]]))
})

test_that("lower-case wells and dashed plate names are parsed", {
  x <- testset.gl
  ids <- indNames(x)
  wells <- paste0(rep(LETTERS[1:8], length.out = length(ids)),
                  rep(1:12, each = 8, length.out = length(ids)))
  capture.output(up <- gl.report.contamination(x, plate = data.frame(id = ids, plate = "1", well = wells),
                                               plot.display = FALSE, verbose = 0))
  capture.output(lo <- gl.report.contamination(x, plate = data.frame(id = ids, plate = "1", well = tolower(wells)),
                                               plot.display = FALSE, verbose = 0))
  expect_equal(lo$pairs, up$pairs)
  x2 <- x
  x2@other$ind.metrics$plate_location <- paste0("PL-1-", wells)
  capture.output(d <- gl.report.contamination(x2, plot.display = FALSE, verbose = 0))
  expect_equal(d$pairs, up$pairs)
  expect_equal(d$ind$well[match(ids, d$ind$id)], wells)
})

test_that("plate wells that do not parse are reported", {
  x <- testset.gl
  plate <- data.frame(id = indNames(x), plate = "1", well = "well")
  out <- capture.output(r <- gl.report.contamination(x, plate = plate, plot.display = FALSE, verbose = 2))
  expect_true(any(grepl("no plate well parses", out)))
  expect_null(r$pairs)
})

# Plate map. Put the dose-pattern host on plate 1 with the rest of population
# A, its partner in the next well (or on plate 2 with the rest of population
# B), and read back the map saved beside the scatter plot.
plate.map <- function(adjacent = TRUE, cols = gl.colors("dis", verbose = 0)) {
  s <- depth.host(function(d) 0.6 * rank(d) / length(d))
  capture.output(r0 <- gl.report.contamination(s$x, plot.display = FALSE,
                                               verbose = 0))
  partner <- r0$ind$top.partner[r0$ind$id == s$host]
  ids <- indNames(s$x)
  on1 <- c(s$host, if (adjacent) partner,
           setdiff(ids[pop(s$x) == "A"], s$host))
  on2 <- setdiff(ids, on1)
  # column-major wells, so the second individual of a plate sits in B1
  wells <- function(n) {
    paste0(LETTERS[(seq_len(n) - 1) %% 8 + 1], (seq_len(n) - 1) %/% 8 + 1)
  }
  plate <- data.frame(id = c(on1, on2),
                      plate = rep(c("1", "2"), c(length(on1), length(on2))),
                      well = c(wells(length(on1)), wells(length(on2))))
  f <- paste0("contamination_map_", as.integer(adjacent))
  capture.output(r <- gl.report.contamination(s$x, plate = plate,
                                              plot.colors.pop = cols,
                                              plot.display = FALSE,
                                              plot.dir = tempdir(),
                                              plot.file = f, verbose = 0))
  c(list(h = r$ind[r$ind$id == s$host, ], partner = partner), read.map(f))
}

# The saved map and the computed data of its rings, arrows and labels
read.map <- function(f) {
  map <- readRDS(file.path(tempdir(), paste0(f, "_plate.RDS")))
  geoms <- vapply(map$layers, function(l) class(l$geom)[1], "")
  ld <- lapply(seq_along(geoms), function(i) ggplot2::layer_data(map, i))
  rings <- Filter(function(d) all(is.na(d$fill)), ld[geoms == "GeomTile"])
  list(map = map,
       arrows = if (any(geoms == "GeomSegment")) ld[[which(geoms == "GeomSegment")]],
       rings = rings[[1]], labels = ld[[which(geoms == "GeomText")]])
}

test_that("the plate map rings a flagged well and draws an arrow from an adjacent partner", {
  m <- plate.map(adjacent = TRUE)
  expect_equal(m$h$flag, "adjacent")
  # only the plate holding a flagged individual is drawn
  expect_equal(as.character(ggplot2::ggplot_build(m$map)$layout$layout$plate),
               "1")
  expect_equal(m$rings$linetype, "solid")
  # one arrow, from the partner's well B1 towards the host's well A1; row A
  # is drawn at the top
  expect_equal(nrow(m$arrows), 1)
  expect_equal(c(m$arrows$x, m$arrows$xend), c(1, 1))
  expect_lt(m$arrows$y, m$arrows$yend)
  # the well is labelled with the host and its foreign-allele percentage
  lab <- m$labels$label
  expect_true(any(grepl(m$h$id, lab)))
  expect_true(any(grepl(paste0(round(100 * m$h$foreign.rate), "%"), lab,
                        fixed = TRUE)))
  # and with its tier 3 pattern, so a hybrid does not read as a contaminant
  expect_equal(m$h$pattern, "dose")
  expect_true(any(grepl("dose", lab, fixed = TRUE)))
})

test_that("the plate map draws a dashed ring and no arrow when the partner is not adjacent", {
  m <- plate.map(adjacent = FALSE)
  expect_equal(m$h$flag, "suspect")
  expect_false(m$h$adjacent)
  expect_equal(m$rings$linetype, "dashed")
  expect_true(is.null(m$arrows) || nrow(m$arrows) == 0)
})

test_that("reciprocal partners in adjacent wells get two arrows side by side", {
  # i1 (A) carries B alleles and i21 (B) carries A alleles at disjoint fixed
  # loci, so each names the other as its partner
  s <- depth.host(function(d) rep(0, length(d)))
  m <- as.matrix(s$x)
  m[1, 1:180] <- 1
  m[21, 181:360] <- 1
  x <- s$x
  x@gen <- new("genlight", m, ploidy = 2)@gen
  indNames(x) <- indNames(s$x); pop(x) <- pop(s$x)
  ids <- indNames(x)
  plate <- data.frame(id = c("i1", "i21", setdiff(ids, c("i1", "i21"))),
                      plate = "1",
                      well = paste0(LETTERS[(0:39) %% 8 + 1], (0:39) %/% 8 + 1))
  f <- "contamination_map_reciprocal"
  capture.output(r <- gl.report.contamination(x, plate = plate,
                                              plot.display = FALSE,
                                              plot.dir = tempdir(),
                                              plot.file = f, verbose = 0))
  expect_equal(r$ind$top.partner[match(c("i1", "i21"), r$ind$id)],
               c("i21", "i1"))
  a <- read.map(f)$arrows
  expect_equal(nrow(a), 2)
  # both arrows run along the A1-B1 edge; they must not share a line
  expect_gt(abs(a$x[1] - a$x[2]), 0.1)
})

test_that("the ring takes the colour of the partner's population", {
  m <- plate.map(adjacent = TRUE)
  # the first layer fills one tile per well of the map data by population
  bg <- ggplot2::layer_data(m$map, 1)
  expect_true(m$partner %in% m$map$data$id)
  fill.b <- unique(bg$fill[m$map$data$pop %in% "B"])
  expect_length(fill.b, 1)
  expect_equal(m$rings$colour, fill.b)
})

test_that("no plate map is saved without plate positions", {
  f <- "contamination_map_none"
  capture.output(gl.report.contamination(testset.gl, plot.display = FALSE,
                                         plot.dir = tempdir(), plot.file = f,
                                         verbose = 0))
  expect_true(file.exists(file.path(tempdir(), paste0(f, ".RDS"))))
  expect_false(file.exists(file.path(tempdir(), paste0(f, "_plate.RDS"))))
})

test_that("plot.colors.pop needs a palette or one colour per population", {
  s <- depth.host(function(d) rep(0.3, length(d)))
  expect_error(capture.output(gl.report.contamination(s$x, plot.colors.pop = "red",
                                                      plot.display = FALSE,
                                                      verbose = 0)),
               "plot.colors.pop")
})

test_that("named plot.colors.pop keeps its names on the plate map", {
  # before: setNames(..., levels(pop(x))) overwrote the names given, so
  # c(B = "red", A = "blue") filled population A red
  m <- plate.map(adjacent = TRUE, cols = c(B = "red", A = "blue"))
  bg <- ggplot2::layer_data(m$map, 1)
  expect_equal(unique(bg$fill[m$map$data$pop %in% "A"]), "blue")
  expect_equal(unique(bg$fill[m$map$data$pop %in% "B"]), "red")
  expect_equal(m$rings$colour, "red")
})

# Read counts for a cohort in which "mix" is 80% the DNA of i2 and 20% that of
# i3 (all three in population A), "mix_R" a second library of the same DNA
# and "i3_R" a second library of i3. Genotypes are called from the reads.
reads.cohort <- function(seed = 1, drop = character(0)) {
  set.seed(seed)
  n.loc <- 3000
  p.a <- runif(n.loc, 0.05, 0.5); p.b <- runif(n.loc, 0.05, 0.5)
  true <- rbind(matrix(rbinom(25 * n.loc, 2, rep(p.a, each = 25)), 25),
                matrix(rbinom(15 * n.loc, 2, rep(p.b, each = 15)), 15))
  rownames(true) <- paste0("i", 1:40)
  # alt-allele fraction of the DNA in each library
  q <- true / 2
  q <- rbind(q, mix = 0.8 * q["i2", ] + 0.2 * q["i3", ])
  q <- rbind(q, mix_R = q["mix", ], i3_R = q["i3", ])
  n <- nrow(q)
  depth <- matrix(rpois(n * n.loc, 40), n)
  # allelic imbalance between libraries, and sequencing error at homozygotes
  qq <- ifelse(q > 0 & q < 1, rbeta(length(q), 30 * q, 30 * (1 - q)),
               pmin(pmax(q, 0.005), 0.995))
  alt <- matrix(rbinom(length(q), depth, qq), n)
  frac <- alt / depth
  g <- ifelse(frac < 0.1, 0, ifelse(frac > 0.9, 2, 1))
  g[depth < 3] <- NA
  keep <- !rownames(q) %in% drop
  x <- new("genlight", g[keep, ], ploidy = 2)
  indNames(x) <- rownames(q)[keep]
  locNames(x) <- paste0("L", seq_len(n.loc))
  pop(x) <- ifelse(rownames(q)[keep] %in% paste0("i", 26:40), "B", "A")
  x <- gl.compliance.check(x, verbose = 0)
  x@other$geno.metrics <- list(ref.count = (depth - alt)[keep, , drop = FALSE],
                               alt.count = alt[keep, , drop = FALSE])
  storage.mode(x@other$geno.metrics$ref.count) <- "integer"
  storage.mode(x@other$geno.metrics$alt.count) <- "integer"
  x
}

test_that("read balance flags a mixed sample and traces both of its sources", {
  x <- reads.cohort()
  capture.output(r <- gl.report.contamination(x, plot.display = FALSE,
                                              verbose = 0))
  ind <- r$ind
  # both libraries of the mixed DNA, and nothing else
  expect_setequal(ind$id[ind$flag == "mixed"], c("mix", "mix_R"))
  h <- ind[ind$id == "mix", ]
  expect_gt(h$skew.z, 3)
  expect_gt(h$read.calls, 1000)
  # the major genome, not the other library of the same DNA
  expect_equal(h$main.source, "i2")
  expect_gt(h$main.kin, 0.2)
  expect_lt(h$extra.share, 0.35)
  # the minor genome; its own replicate is not counted as the runner-up
  expect_true(h$second.source %in% c("i3", "i3_R"))
  expect_false(h$second.next %in% c("i3", "i3_R"))
  expect_gt(h$second.gap, 0.15)
  # 20% of the DNA came from i3
  expect_gt(h$mix.fraction, 0.12)
  expect_lt(h$mix.fraction, 0.3)
  # the sources themselves read as balanced
  expect_true(all(ind$skew.z[ind$id %in% c("i2", "i3", "i3_R")] < 3))
  expect_true(all(is.na(ind$main.source[ind$flag != "mixed"])))
})

test_that("read balance does not trace sources when the major genome is absent", {
  # without i2 the most related individual is not the sample's major genome,
  # so the extra sites hold the sample's own heterozygotes
  x <- reads.cohort(drop = c("i2", "mix_R"))
  capture.output(r <- gl.report.contamination(x, plot.display = FALSE,
                                              verbose = 0))
  h <- r$ind[r$ind$id == "mix", ]
  expect_equal(h$flag, "mixed")
  expect_gte(h$extra.share, 0.35)
  expect_true(is.na(h$second.source) && is.na(h$mix.fraction))
})

test_that("read balance skips individuals with few deep genotypes", {
  x <- reads.cohort()
  # leave mix with 500 genotypes of 20 or more reads
  i <- match("mix", indNames(x))
  low <- seq_len(nLoc(x)) > 500
  x@other$geno.metrics$ref.count[i, low] <- 5L
  x@other$geno.metrics$alt.count[i, low] <- 5L
  capture.output(r <- gl.report.contamination(x, plot.display = FALSE,
                                              verbose = 0))
  h <- r$ind[r$ind$id == "mix", ]
  expect_lte(h$read.calls, 500)
  expect_true(is.na(h$skew.z))
  expect_false(h$flag == "mixed")
})

test_that("without read counts the read-balance columns are NA", {
  x <- reads.cohort()
  x@other$geno.metrics <- NULL
  out <- capture.output(r <- gl.report.contamination(x, plot.display = FALSE,
                                                     verbose = 2))
  expect_true(any(grepl("read.counts = TRUE", out)))
  expect_true(all(is.na(r$ind$skew.z)))
  expect_false(any(r$ind$flag == "mixed"))
  expect_true(all(c("read.calls", "skew.rate", "skew.z", "main.source",
                    "main.kin", "n.extra", "extra.share", "second.source",
                    "second.share", "second.next", "second.gap",
                    "mix.fraction") %in% names(r$ind)))
})

test_that("gl.filter.contamination removes mixed samples by default", {
  x <- reads.cohort()
  capture.output(y <- gl.filter.contamination(x, verbose = 0))
  expect_false(any(c("mix", "mix_R") %in% indNames(y)))
  capture.output(y2 <- gl.filter.contamination(x, flag = "suspect",
                                               verbose = 0))
  expect_true(all(c("mix", "mix_R") %in% indNames(y2)))
})
