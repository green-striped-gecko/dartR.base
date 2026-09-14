# Characterization tests for gl.report.shannon
# Baseline snapshotted before review (review-gl.report.shannon), then
# updated for the approved findings of the function review:
#   [approved F1] level = "beta"/"gamma" now compute the real population-level
#     partition (Ma, Li & Zhang 2020 matrix forms) instead of the degenerate
#     beta == 1 / gamma == alpha per-individual constants.
#   [approved F2] verbose = NULL default: the global gl.set.verbosity setting
#     is honoured.
#   [approved F3] verbose = 0 is fully silent (melt message gone, plot gated).
#   [approved F4] level and order are validated with informative errors.
#   [approved F6] individuals with no non-missing, non-zero dosages return NA
#     rows (was 0/1/Inf) with a gated warning.
# The alpha profile itself is pinned unchanged: those values predate the
# review and were verified against independent hand computation.

# Small local fixture: 3 individuals x 4 loci, dosages 0/1/2, one
# all-reference individual (i1).
make_small <- function() {
  m <- matrix(c(0, 0, 0, 0,
                1, 2, 0, 1,
                2, 2, 1, 0), nrow = 3, byrow = TRUE)
  rownames(m) <- paste0("i", 1:3)
  colnames(m) <- paste0("L", 1:4)
  g <- new("genlight", m)
  pop(g) <- factor(rep("a", 3))
  ploidy(g) <- rep(2L, 3)
  g
}

# Two-population fixture: 4 informative individuals x 4 loci, pops A and B.
make_twopop <- function() {
  m <- matrix(c(1, 2, 0, 1,
                2, 2, 1, 0,
                0, 1, 1, 2,
                2, 0, 2, 1), nrow = 4, byrow = TRUE)
  rownames(m) <- paste0("i", 1:4)
  colnames(m) <- paste0("L", 1:4)
  g <- new("genlight", m)
  pop(g) <- factor(c("A", "A", "B", "B"))
  ploidy(g) <- rep(2L, 4)
  g
}

# Independent Hill-number implementations (hand computation, no d.chao):
# vector form for a single abundance vector, matrix forms for the Ma, Li &
# Zhang (2020) population partition (rows = individuals, cols = loci).
hill_vec <- function(v, q) {
  v <- v[!is.na(v)]; v <- v[v > 0]; p <- v / sum(v)
  if (q == 1) exp(-sum(p * log(p))) else (sum(p^q))^(1 / (1 - q))
}
hill_gamma <- function(A, q) hill_vec(colSums(A), q)      # pooled abundances
hill_alpha <- function(A, q) {                            # normalised mean
  eA <- A / sum(A); eA <- eA[eA > 0]; N <- nrow(A)
  if (q == 1) exp(-sum(eA * log(eA)) - log(N))
  else (1 / N) * (sum(eA^q))^(1 / (1 - q))
}

quiet_shannon <- function(...) {
  res <- NULL
  suppressMessages(
    invisible(capture.output(
      res <- gl.report.shannon(..., plot.display = FALSE, verbose = 0)
    ))
  )
  res
}

num <- function(r) {
  n <- as.matrix(r[, -1]); mode(n) <- "numeric"; unname(n)
}

test_that("alpha profile on testset.gl: structure and pinned values", {
  r <- quiet_shannon(testset.gl)
  expect_s3_class(r, "data.frame")
  expect_equal(dim(r), c(250, 6))                       # 250 ind, ID + q0..q4
  expect_equal(colnames(r), c("ID", paste0("q", 0:4)))
  expect_equal(r$ID, indNames(testset.gl))
  n <- num(r)
  expect_equal(n[1, ],                                  # AA010915
               c(70, 69.4005218514, 69.0149253731, 68.7626406698,
                 68.5935462245), tolerance = 1e-9)
  expect_equal(n[250, ],                                # AA001450
               c(84, 83.3980790444, 83.0123456790, 82.7604490410,
                 82.5917730496), tolerance = 1e-9)
  expect_equal(unname(colMeans(n)),
               c(76.3960000000, 75.8498421158, 75.4982703728,
                 75.2681132211, 75.1137969530), tolerance = 1e-9)
})

test_that("alpha values equal hand-computed Hill numbers on dosages", {
  # Independent recomputation: per individual, drop NA and zero dosages,
  # p = dosage/sum(dosage), qD = (sum p^q)^(1/(1-q)), q1 = exp(Shannon).
  r <- quiet_shannon(testset.gl)
  mat <- as.matrix(testset.gl)
  hand <- t(sapply(seq_len(nrow(mat)), function(i)
    sapply(0:4, function(q) hill_vec(mat[i, ], q))))
  expect_equal(num(r), unname(hand), tolerance = 1e-12)
})

test_that("beta and gamma report the real population partition [approved F1]", {
  g <- make_twopop()
  rg <- quiet_shannon(g, level = "gamma")
  rb <- quiet_shannon(g, level = "beta")
  expect_equal(colnames(rg), c("pop", paste0("q", 0:4)))
  expect_equal(rg$pop, c("A", "B"))
  ng <- num(rg); nb <- num(rb)
  # Hand-verified gamma for pop A: pooled abundances (3, 4, 1, 1), tot 9.
  pA <- c(3, 4, 1, 1) / 9
  expect_equal(ng[1, 1], 4)                              # pooled richness
  expect_equal(ng[1, 2], exp(-sum(pA * log(pA))), tolerance = 1e-12)
  expect_equal(ng[1, 3], 1 / sum(pA^2), tolerance = 1e-12)
  # Gamma independently from pooled frequencies, both populations, all orders
  mat <- as.matrix(g); mat[is.na(mat)] <- 0
  for (i in 1:2) {
    A <- mat[pop(g) == c("A", "B")[i], , drop = FALSE]
    expect_equal(ng[i, ], sapply(0:4, function(q) hill_gamma(A, q)),
                 tolerance = 1e-12)
    # Partition identity gamma = alpha x beta at every order, alpha computed
    # independently of the function
    expect_equal(ng[i, ], sapply(0:4, function(q) hill_alpha(A, q)) * nb[i, ],
                 tolerance = 1e-12)
  }
  # Beta no longer the degenerate constant 1; bounded by (1, N]
  expect_true(all(nb > 1) && all(nb <= 2))
})

test_that("partition identity gamma = alpha x beta on testset.gl [approved F1]", {
  rg <- quiet_shannon(testset.gl, level = "gamma")
  rb <- quiet_shannon(testset.gl, level = "beta")
  expect_equal(rg$pop, levels(pop(testset.gl)))
  ng <- num(rg); nb <- num(rb)
  mat <- as.matrix(testset.gl); mat[is.na(mat)] <- 0
  informative <- rowSums(mat) > 0
  for (i in seq_along(levels(pop(testset.gl)))) {
    A <- mat[pop(testset.gl) == levels(pop(testset.gl))[i] & informative, ,
             drop = FALSE]
    alpha_i <- sapply(0:4, function(q) hill_alpha(A, q))
    expect_equal(ng[i, ], sapply(0:4, function(q) hill_gamma(A, q)),
                 tolerance = 1e-12)
    expect_equal(ng[i, ], alpha_i * nb[i, ], tolerance = 1e-12)
  }
  expect_false(all(nb == 1))                             # was the F1 defect
})

test_that("report contract: input untouched, invisible return, plot decoupled", {
  before <- serialize(testset.gl, NULL)
  suppressMessages(invisible(capture.output(
    v <- withVisible(gl.report.shannon(testset.gl, plot.display = FALSE,
                                       verbose = 0))
  )))
  expect_identical(serialize(testset.gl, NULL), before)
  expect_false(v$visible)
  r_off <- v$value
  pdf(NULL)
  suppressMessages(invisible(capture.output(
    r_on <- gl.report.shannon(testset.gl, plot.display = TRUE, verbose = 2)
  )))
  dev.off()
  expect_identical(r_on, r_off)                          # PLT3
})

test_that("alpha results independent of population assignment", {
  g <- testset.gl
  set.seed(42)
  pop(g) <- factor(sample(letters[1:5], nInd(g), replace = TRUE))
  expect_identical(quiet_shannon(g), quiet_shannon(testset.gl))
})

test_that("verbose = 0 is fully silent, stdout and messages [approved F3]", {
  msgs <- character(0)
  o <- capture.output(
    withCallingHandlers(
      r <- gl.report.shannon(testset.gl, plot.display = FALSE, verbose = 0),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
      }
    ),
    type = "output"
  )
  expect_length(o, 0)
  expect_length(msgs, 0)
})

test_that("global verbosity 0 is honoured (verbose = NULL) [approved F2]", {
  op <- options(dartR_verbose = 0); on.exit(options(op))
  o <- capture.output(
    suppressMessages(r <- gl.report.shannon(testset.gl, plot.display = FALSE)),
    type = "output"
  )
  expect_length(o, 0)
})

test_that("all-reference individual: NA row and gated warning [approved F6]", {
  g <- make_small()
  r <- quiet_shannon(g)
  n <- num(r)
  expect_true(all(is.na(n[1, ])))                        # was 0/1/Inf
  expect_equal(n[2, 1], 3)                               # i2: three non-zero loci
  o <- capture.output(
    r1 <- gl.report.shannon(g, plot.display = FALSE, verbose = 1),
    type = "output"
  )
  expect_true(any(grepl("i1", o)))                       # warning names it
  o0 <- capture.output(
    r0 <- gl.report.shannon(g, plot.display = FALSE, verbose = 0),
    type = "output"
  )
  expect_length(o0, 0)                                   # gated at verbose 0
})

test_that("degenerate individuals are excluded from the partition [approved F1, F6]", {
  g <- make_small()                                      # i1 all-reference
  rg <- quiet_shannon(g, level = "gamma")
  mat <- as.matrix(g)[2:3, ]; mat[is.na(mat)] <- 0
  expect_equal(num(rg)[1, ], sapply(0:4, function(q) hill_gamma(mat, q)),
               tolerance = 1e-12)
})

test_that("order = 1 returns a single q0 column; single individual works", {
  g <- make_small()
  r1 <- quiet_shannon(g[2:3, ], order = 1)
  expect_equal(colnames(r1), c("ID", "q0"))
  rs <- quiet_shannon(g[2, , drop = FALSE])
  expect_equal(nrow(rs), 1)
})

test_that("invalid level and order fail fast with informative errors [approved F4]", {
  g <- make_small()
  expect_error(quiet_shannon(g[2:3, ], level = "banana"),
               "level must be one of")
  expect_error(quiet_shannon(g[2:3, ], order = 0),
               "order must be a positive whole number")
  expect_error(quiet_shannon(g[2:3, ], order = 2.5),
               "order must be a positive whole number")
})

test_that("SilicoDArT data are rejected (accept = 'SNP')", {
  expect_error(quiet_shannon(testset.gs),
               "found SilicoDArT expecting SNP")
})
