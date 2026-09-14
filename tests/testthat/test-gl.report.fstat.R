# Characterization tests for gl.report.fstat
# Baseline snapshotted before review (review-gl.report.fstat), on
# dartR.data 1.2.5, then updated when the approved findings were applied.
# Assertions that changed carry an "[approved Fn]" comment naming the
# finding that changed them; every other assertion is untouched.
#
# Point estimates are pinned against utils.basic.stats rather than against
# literals wherever the literal depends on which version of that helper is
# in the tree: PR #309 (the utils.basic.stats review) is open but not yet
# merged into dev, and it moves Fst and Gst_H by 2e-4 on platypus.gl. The
# gl.report.fstat review does not touch the point-estimate path.

# ---------------------------------------------------------------- helpers

# replace the genotype matrix of a compliant object, keeping metadata
mk_fixture <- function(template, m) {
  g <- template[, seq_len(ncol(m))]
  g@gen <- new("genlight", m, ploidy = 2)@gen
  g@n.loc <- ncol(m)
  g
}

two_pop <- function() {
  x <- platypus.gl[pop(platypus.gl) %in%
                     c("SEVERN_ABOVE", "SEVERN_BELOW"), ]
  pop(x) <- droplevels(pop(x))
  x
}

# ------------------------------------------------- report contract (FS8)

test_that("input object comes back untouched and no history is appended", {
  xcopy <- platypus.gl
  nh <- length(platypus.gl@other$history)
  invisible(gl.report.fstat(platypus.gl, verbose = 0, plot.display = FALSE))
  expect_identical(xcopy, platypus.gl)
  expect_equal(length(platypus.gl@other$history), nh)
})

test_that("results are independent of plotting (PLT3)", {
  a <- gl.report.fstat(platypus.gl, verbose = 0, plot.display = TRUE)
  b <- gl.report.fstat(platypus.gl, verbose = 0, plot.display = FALSE)
  expect_identical(a, b)
})

test_that("verbose = 0 is fully silent and draws nothing", {
  # [approved F4] plot.display is forced FALSE at verbose 0, so the heatmap
  # path is never entered. Before the fix this emitted one line: the
  # dependency guard printed from inside gl.plot.heatmap. The absence of
  # that line is also the evidence that no graphic is attempted.
  expect_length(
    capture.output(invisible(
      gl.report.fstat(platypus.gl, verbose = 0, plot.display = TRUE))),
    0)
  expect_length(
    capture.output(invisible(
      gl.report.fstat(platypus.gl, verbose = 0, plot.display = FALSE))),
    0)
})

test_that("verbose 3 prints the results summary and verbose 2 does not", {
  # [approved F9] the summary moved from verbose >= 2 to verbose >= 3, so
  # levels 2 and 3 are no longer indistinguishable.
  o1 <- capture.output(invisible(
    gl.report.fstat(platypus.gl, verbose = 1, plot.display = FALSE)))
  o2 <- capture.output(invisible(
    gl.report.fstat(platypus.gl, verbose = 2, plot.display = FALSE)))
  o3 <- capture.output(invisible(
    gl.report.fstat(platypus.gl, verbose = 3, plot.display = FALSE)))
  expect_false(identical(o2, o3))
  expect_lt(length(o2), 10)
  expect_gt(length(o3), 40)
  expect_lte(length(o1), length(o2))
})

# ------------------------------------------------------ pinned statistics

test_that("platypus.gl pairwise statistics are pinned", {
  r <- gl.report.fstat(platypus.gl, verbose = 0, plot.display = FALSE)
  # Fstp and Dest are identical under both versions of utils.basic.stats
  expect_equal(unname(r[[1]]$Fstp[lower.tri(r[[1]]$Fstp)]),
               c(0.0599, 0.0754, 0.0836))
  expect_equal(unname(r[[1]]$Dest[lower.tri(r[[1]]$Dest)]),
               c(0.0101, 0.0134, 0.0150))
  # Fst and Gst_H are 0.0309/0.0392/0.0436 and 0.0753/0.0958/0.1062 with
  # dev's utils.basic.stats, and 0.0307/0.0392/0.0434 and
  # 0.0754/0.0958/0.1063 with the reviewed one (PR #309). Pin the
  # invariant instead: every reported cell is exactly the helper's value.
  pops <- seppop(platypus.gl)
  prs <- t(combn(length(pops), 2))
  for (k in seq_len(nrow(prs))) {
    tp <- rbind.dartR(pops[[prs[k, 1]]], pops[[prs[k, 2]]])
    d <- utils.basic.stats(tp)$overall[c("Fst", "Fstp", "Dest", "Gst_H")]
    for (s in c("Fst", "Fstp", "Dest", "Gst_H")) {
      expect_equal(r[[1]][[s]][prs[k, 2], prs[k, 1]], unname(d[[s]]))
    }
  }
})

test_that("point estimates match hierfstat::basic.stats (Nei 1987)", {
  skip_if_not_installed("hierfstat")
  r <- gl.report.fstat(platypus.gl, verbose = 0, plot.display = FALSE)
  pops <- seppop(platypus.gl)
  prs <- t(combn(length(pops), 2))
  for (k in seq_len(nrow(prs))) {
    tp <- rbind.dartR(pops[[prs[k, 1]]], pops[[prs[k, 2]]])
    hf <- suppressWarnings(
      hierfstat::basic.stats(
        hierfstat::genind2hierfstat(gl2gi(tp, verbose = 0))))
    ov <- hf$overall
    # Fstp and Dest agree exactly at the function's own 4-decimal rounding
    expect_equal(r[[1]]$Fstp[prs[k, 2], prs[k, 1]],
                 round(unname(ov["Fstp"]), 4), tolerance = 1e-4)
    expect_equal(r[[1]]$Dest[prs[k, 2], prs[k, 1]],
                 round(unname(ov["Dest"]), 4), tolerance = 1e-4)
    # Fst is exact once PR #309 lands; with dev's utils.basic.stats two of
    # the three pairs sit 2e-4 away, from the unguarded harmonic mean of
    # per-population sample sizes
    expect_lt(abs(r[[1]]$Fst[prs[k, 2], prs[k, 1]] -
                    round(unname(ov["Fst"]), 4)), 3e-4)
  }
})

test_that("the reported statistic is Nei's Gst, not Weir & Cockerham", {
  skip_if_not_installed("hierfstat")
  pops <- seppop(platypus.gl)
  tp <- rbind.dartR(pops[[1]], pops[[2]])
  wc <- hierfstat::wc(
    suppressWarnings(hierfstat::genind2hierfstat(gl2gi(tp, verbose = 0))))
  r <- gl.report.fstat(platypus.gl, verbose = 0, plot.display = FALSE)
  # Fst (biased Nei Gst) is materially below theta
  expect_gt(abs(r[[1]]$Fst[2, 1] - wc$FST), 0.02)
  # Fstp (Nei's unbiased Gst) sits within 0.001 of theta on this data
  expect_lt(abs(r[[1]]$Fstp[2, 1] - wc$FST), 0.001)
})

test_that("label integrity holds for non-alphabetical population levels", {
  x2 <- platypus.gl
  pop(x2) <- factor(as.character(pop(platypus.gl)),
                    levels = c("TENTERFIELD", "SEVERN_BELOW",
                               "SEVERN_ABOVE"))
  r1 <- gl.report.fstat(platypus.gl, verbose = 0, plot.display = FALSE)
  r2 <- gl.report.fstat(x2, verbose = 0, plot.display = FALSE)
  nm <- c("SEVERN_ABOVE", "SEVERN_BELOW", "TENTERFIELD")
  expect_equal(r1[[1]]$Fst[nm, nm], r2[[1]]$Fst[nm, nm])
  expect_equal(r1[[1]]$Fst, t(r1[[1]]$Fst))
})

test_that("matrix cells agree with a direct per-pair computation", {
  r <- gl.report.fstat(possums.gl, verbose = 0, plot.display = FALSE)
  pp <- seppop(possums.gl)
  for (pr in list(c(1, 4), c(2, 7), c(3, 10))) {
    tp <- rbind.dartR(pp[[pr[1]]], pp[[pr[2]]])
    direct <- unname(utils.basic.stats(tp)$overall["Fstp"])
    expect_equal(unname(r[[1]]$Fstp[names(pp)[pr[1]], names(pp)[pr[2]]]),
                 direct)
  }
  expect_equal(r[[1]]$Fstp["A", "B"], 0.2757)
  expect_equal(r[[1]]$Fstp["A", "J"], 0.2917)
})

# ------------------------------------------------------- return structure

test_that("return shape varies by nboots and number of populations", {
  # > 2 populations, no bootstrap: list(Stat_matrices, Stat_tables)
  r <- gl.report.fstat(possums.gl, verbose = 0, plot.display = FALSE)
  expect_type(r, "list")
  expect_length(r, 2)
  # [approved F6] the second element is named
  expect_identical(names(r), c("Stat_matrices", "Stat_tables"))
  expect_identical(names(r[[1]]), c("Fst", "Fstp", "Dest", "Gst_H"))
  expect_equal(dim(r[[1]]$Fst), c(10L, 10L))
  # [approved F11] columns are the bare pair names in every branch
  expect_false(any(grepl("^Stat_tables\\.", colnames(r[[2]]))))
  expect_identical(colnames(r[[2]])[1], "A_vs_B")
  expect_equal(dim(r[[2]]), c(4L, 45L))

  # 2 populations, no bootstrap: a bare data.frame, not "two lists"
  r2 <- gl.report.fstat(two_pop(), verbose = 0, plot.display = FALSE)
  expect_s3_class(r2, "data.frame")
  expect_equal(dim(r2), c(4L, 1L))
  expect_identical(rownames(r2), c("Fst", "Fstp", "Dest", "Gst_H"))
  expect_identical(names(r2), "SEVERN_ABOVE_vs_SEVERN_BELOW")

  # with bootstrap
  set.seed(2)
  rb <- gl.report.fstat(two_pop(), nboots = 50, CI.type = "perc",
                        verbose = 0, plot.display = FALSE)
  expect_identical(names(rb), c("Stat_tables", "Confidence_Intervals"))
  expect_identical(colnames(rb$Confidence_Intervals),
                   c("Value", "LCI", "HCI"))

  set.seed(2)
  rb3 <- gl.report.fstat(platypus.gl, nboots = 50, CI.type = "perc",
                         verbose = 0, plot.display = FALSE)
  expect_identical(names(rb3), c("Stat_matrices", "Confidence_Intervals"))
  expect_length(rb3$Confidence_Intervals, 3)
  expect_identical(names(rb3$Confidence_Intervals),
                   c("SEVERN_ABOVE_vs_SEVERN_BELOW",
                     "SEVERN_ABOVE_vs_TENTERFIELD",
                     "SEVERN_BELOW_vs_TENTERFIELD"))
})

# ------------------------------------------------------------- bootstrap

test_that("boot() draws loci, and every locus can be drawn", {
  # [approved F1, BLOCKER] boot::boot is handed a one-column frame of locus
  # positions, so its row indices index loci. Before the fix it was handed
  # an individuals-by-loci frame and the statistic applied the row indices
  # to the columns, capping every replicate at the first nInd loci.
  pops <- seppop(platypus.gl)
  tp <- rbind.dartR(pops[[1]], pops[[2]])
  seen <- new.env(); seen$v <- integer(0)
  spy <- function(loc.index, indices, gen.mat, pops.info) {
    seen$v <- c(seen$v, loc.index$loc[indices]); c(0, 0, 0, 0)
  }
  set.seed(1)
  invisible(boot::boot(data = data.frame(loc = seq_len(nLoc(tp))),
                       statistic = spy, gen.mat = as.matrix(tp),
                       pops.info = as.character(pop(tp)), R = 20))
  expect_equal(nInd(tp), 40L)
  expect_equal(nLoc(tp), 1000L)
  # replicates are nLoc long, not nInd long, and reach the last locus
  expect_equal(length(seen$v) %% nLoc(tp), 0)
  expect_equal(max(seen$v), nLoc(tp))
  expect_equal(length(unique(seen$v)), nLoc(tp))
})

test_that("bootstrap CIs bracket the reported point estimate", {
  # [approved F1, BLOCKER] fixture: the first nInd loci are undifferentiated,
  # the remaining loci are fixed differences. Before the fix every replicate
  # used only the first nInd loci, giving Value 0.8622 with a 95 per cent
  # percentile interval of [-0.0175, -0.0175].
  x <- platypus.gl[pop(platypus.gl) %in%
                     c("SEVERN_ABOVE", "SEVERN_BELOW"), 1:300]
  pop(x) <- droplevels(pop(x))
  m <- as.matrix(x)
  ni <- nrow(m)
  grp <- as.character(pop(x))
  for (j in seq_len(ni)) m[, j] <- rep(c(0L, 1L, 2L), length.out = ni)
  for (j in (ni + 1):ncol(m)) m[, j] <- ifelse(grp == grp[1], 0L, 2L)
  fx <- mk_fixture(x, m)
  set.seed(7)
  r <- gl.report.fstat(fx, nboots = 200, CI.type = "perc", verbose = 0,
                       plot.display = FALSE)
  ci <- r$Confidence_Intervals
  expect_equal(ci[["Value"]][1], 0.8622)
  expect_gt(ci[["LCI"]][1], 0.75)
  expect_true(all(ci[["LCI"]] <= ci[["Value"]] &
                    ci[["Value"]] <= ci[["HCI"]]))
  # the replicates vary; they are not a single degenerate value
  expect_gt(ci[["HCI"]][1] - ci[["LCI"]][1], 0.01)
})

test_that("the bootstrap and the point estimate share one estimator", {
  # [approved F2, HIGH] the inlined copy of utils.basic.stats was deleted;
  # the statistic now calls the same helper the point estimate calls. A
  # replicate built from the identity index therefore reproduces the point
  # estimate exactly.
  pops <- seppop(platypus.gl)
  tp <- rbind.dartR(pops[[1]], pops[[2]])
  gm <- as.matrix(tp)
  identity_rep <- utils.basic.stats(
    new("genlight", gen = gm[, seq_len(ncol(gm)), drop = FALSE],
        ploidy = 2, pop = as.character(pop(tp)),
        parallel = FALSE))$overall[c("Fst", "Fstp", "Dest", "Gst_H")]
  point <- utils.basic.stats(tp)$overall[c("Fst", "Fstp", "Dest", "Gst_H")]
  expect_identical(unname(identity_rep), unname(point))

  # [approved F2] loci absent from one population of a pair used to abort
  # the whole call from inside boot.ci; the shared estimator completes.
  y <- platypus.gl[pop(platypus.gl) %in%
                     c("SEVERN_ABOVE", "SEVERN_BELOW"), 1:200]
  pop(y) <- droplevels(pop(y))
  m <- as.matrix(y)
  m[which(as.character(pop(y)) == "SEVERN_ABOVE"), 1:40] <- NA
  fx <- mk_fixture(y, m)
  cur <- utils.basic.stats(fx)$overall[c("Fst", "Fstp", "Dest", "Gst_H")]
  r0 <- gl.report.fstat(fx, verbose = 0, plot.display = FALSE)
  expect_equal(unname(unlist(r0[, 1])), unname(cur))
  set.seed(3)
  rb <- gl.report.fstat(fx, nboots = 200, CI.type = "perc", verbose = 0,
                        plot.display = FALSE)
  expect_equal(unname(rb$Confidence_Intervals$Value), unname(cur))
  expect_false(any(is.na(rb$Confidence_Intervals$LCI)))
})

test_that("bootstrap is reproducible only via the global RNG", {
  # The function takes no seed argument; the roxygen states the set.seed
  # guidance. The boot.ci "extreme order statistics" warnings raised at
  # nboots = 30 are a consequence of 30 being too few replicates for a
  # 95 per cent percentile interval, not of the F1 defect: a correct locus
  # bootstrap raises them too, and they disappear at 200 replicates.
  set.seed(11)
  a <- suppressWarnings(
    gl.report.fstat(platypus.gl, nboots = 30, CI.type = "perc",
                    verbose = 0, plot.display = FALSE))
  set.seed(11)
  b <- suppressWarnings(
    gl.report.fstat(platypus.gl, nboots = 30, CI.type = "perc",
                    verbose = 0, plot.display = FALSE))
  expect_identical(a, b)
  set.seed(11)
  expect_no_warning(
    gl.report.fstat(platypus.gl, nboots = 200, CI.type = "perc",
                    verbose = 0, plot.display = FALSE))
})

test_that("conf and CI.type are honoured", {
  set.seed(5)
  w1 <- gl.report.fstat(platypus.gl, nboots = 200, conf = 0.50,
                        CI.type = "perc", verbose = 0, plot.display = FALSE)
  set.seed(5)
  w2 <- gl.report.fstat(platypus.gl, nboots = 200, conf = 0.99,
                        CI.type = "perc", verbose = 0, plot.display = FALSE)
  d1 <- diff(as.numeric(w1$Confidence_Intervals[[1]][1, 2:3]))
  d2 <- diff(as.numeric(w2$Confidence_Intervals[[1]][1, 2:3]))
  expect_gt(d2, d1)
  set.seed(3)
  n1 <- gl.report.fstat(platypus.gl, nboots = 60, CI.type = "norm",
                        verbose = 0, plot.display = FALSE)
  set.seed(3)
  p1 <- gl.report.fstat(platypus.gl, nboots = 60, CI.type = "perc",
                        verbose = 0, plot.display = FALSE)
  expect_false(isTRUE(all.equal(n1$Confidence_Intervals[[1]][1, 2],
                                p1$Confidence_Intervals[[1]][1, 2])))
})

test_that("bad CI.type and bad nboots fail with informative messages", {
  # [approved F5] each argument is validated before any work is done
  expect_error(
    gl.report.fstat(platypus.gl, nboots = 200, CI.type = "rubbish",
                    verbose = 0, plot.display = FALSE),
    "CI.type must be one of")
  expect_error(
    gl.report.fstat(platypus.gl, nboots = 1, verbose = 0,
                    plot.display = FALSE),
    "nboots = 1 gives a bootstrap distribution of a single value")
  expect_error(
    gl.report.fstat(platypus.gl, nboots = -5, CI.type = "perc",
                    verbose = 0, plot.display = FALSE),
    "nboots must be a single non-negative whole number")
  expect_error(
    gl.report.fstat(platypus.gl, nboots = 0.5, CI.type = "perc",
                    verbose = 0, plot.display = FALSE),
    "nboots must be a single non-negative whole number")
  expect_error(
    gl.report.fstat(platypus.gl, conf = 95, verbose = 0,
                    plot.display = FALSE),
    "conf must be a single number")
})

test_that("CI.type = 'bca' is refused below the documented minimum", {
  # [approved F8] 30 replicates used to abort inside boot.ci with
  # "estimated adjustment 'a' is NA"; the default CI.type is unchanged.
  expect_error(
    gl.report.fstat(platypus.gl, nboots = 30, verbose = 0,
                    plot.display = FALSE),
    "requires at least 200 bootstrap replicates")
})

# ------------------------------------------------------------ dispatch

test_that("SilicoDArT data is refused", {
  # [approved F3, DAT7] presence/absence scores were read as SNP dosages:
  # a score of 1 was counted as a heterozygote and frequencies halved.
  expect_error(
    gl.report.fstat(testset.gs, verbose = 0, plot.display = FALSE),
    "found SilicoDArT expecting SNP")
})

test_that("FBM-backed objects are handled", {
  skip_if_not(exists("gl.gen2fbm"))
  fb <- gl.gen2fbm(platypus.gl, verbose = 0)
  a <- gl.report.fstat(fb, verbose = 0, plot.display = FALSE)
  b <- gl.report.fstat(platypus.gl, verbose = 0, plot.display = FALSE)
  expect_equal(a[[1]], b[[1]])
})

# ---------------------------------------------------------- edge cases

test_that("populations of one individual are dropped and announced", {
  # [approved F10, VRB4] the exclusion changes the result, so it is
  # announced from verbose 1 and the dropped populations are named.
  si <- platypus.gl
  p <- as.character(pop(si)); p[1] <- "SOLO"; pop(si) <- factor(p)
  o0 <- capture.output(invisible(
    gl.report.fstat(si, verbose = 0, plot.display = FALSE)))
  o1 <- capture.output(invisible(
    gl.report.fstat(si, verbose = 1, plot.display = FALSE)))
  o2 <- capture.output(invisible(
    gl.report.fstat(si, verbose = 2, plot.display = FALSE)))
  expect_false(any(grepl("more than one", o0)))
  expect_true(any(grepl("more than one", o1)))
  expect_true(any(grepl("SOLO", o1)))
  expect_true(any(grepl("more than one", o2)))
  r <- gl.report.fstat(si, verbose = 0, plot.display = FALSE)
  expect_equal(nrow(r[[1]]$Fst), 3L)
})

test_that("a single population fails with an informative message", {
  # [approved F5] checked precondition rather than combn(1, 2)
  one <- platypus.gl[pop(platypus.gl) == "SEVERN_ABOVE", ]
  pop(one) <- droplevels(pop(one))
  expect_error(gl.report.fstat(one, verbose = 0, plot.display = FALSE),
               "at least two populations of more than one individual")
})

test_that("data of only singleton populations fails informatively", {
  # [approved F5] previously surfaced as gl.keep.pop's internal message
  allsi <- platypus.gl[1:4, ]
  pop(allsi) <- factor(paste0("P", 1:4))
  expect_error(gl.report.fstat(allsi, verbose = 0, plot.display = FALSE),
               "at least two populations of more than one individual")
})

test_that("monomorphic-only data returns NaN rather than erroring", {
  m <- as.matrix(platypus.gl[, 1:50]); m[, ] <- 0L
  mono <- mk_fixture(platypus.gl, m)
  r <- gl.report.fstat(mono, verbose = 0, plot.display = FALSE)
  expect_true(all(is.nan(r[[1]]$Fst[lower.tri(r[[1]]$Fst)])))
})

test_that("populations with identical frequencies give small negative Fst", {
  m <- as.matrix(platypus.gl[, 1:50])
  m[, ] <- rep(c(0L, 1L, 2L), length.out = length(m))
  idn <- mk_fixture(platypus.gl, m)
  r <- gl.report.fstat(idn, verbose = 0, plot.display = FALSE)
  v <- r[[1]]$Fst[lower.tri(r[[1]]$Fst)]
  expect_true(all(v < 0))
  expect_true(all(abs(v) < 0.02))
})

test_that("an invalid plot.stat is refused before any work is done", {
  # [approved F5, F7] plot.stat used to be dereferenced lazily inside the
  # heatmap, so a typo surfaced as "your 'Plots' pane is too small".
  expect_error(
    gl.report.fstat(platypus.gl, verbose = 0, plot.stat = "Nonsense",
                    plot.display = TRUE),
    "plot.stat must be one of")
})

# ------------------------------------------- cross-implementation pins

test_that("Fstp tracks StAMPP's Weir & Cockerham Fst used by gl.fst.pop", {
  skip_if_not_installed("StAMPP")
  a <- gl.report.fstat(possums.gl, verbose = 0, plot.display = FALSE)
  b <- as.matrix(as.dist(gl.fst.pop(possums.gl, nboots = 1, verbose = 0)))
  nm <- rownames(a[[1]]$Fstp)
  b <- b[nm, nm]
  ut <- upper.tri(b)
  # Nei's unbiased Gst and W&C theta agree to 1e-4 on possums.gl
  expect_lt(max(abs(a[[1]]$Fstp[ut] - b[ut])), 1e-4)
  # the biased Nei Gst does not
  expect_gt(max(abs(a[[1]]$Fst[ut] - b[ut])), 0.1)
})

test_that("Fstp and gl.fst.pop stay close under missing data", {
  skip_if_not_installed("StAMPP")
  m <- as.matrix(platypus.gl[, 1:300])
  set.seed(4)
  m[sample(length(m), floor(0.30 * length(m)))] <- NA
  m[which(as.character(pop(platypus.gl)) == "SEVERN_ABOVE"), 1:30] <- NA
  fx <- mk_fixture(platypus.gl, m)
  a <- gl.report.fstat(fx, verbose = 0, plot.display = FALSE)
  b <- as.matrix(as.dist(gl.fst.pop(fx, nboots = 1, verbose = 0)))
  nm <- rownames(a[[1]]$Fstp); b <- b[nm, nm]; ut <- upper.tri(b)
  expect_false(any(is.na(a[[1]]$Fstp[ut])))
  expect_false(any(is.na(b[ut])))
  expect_lt(max(abs(a[[1]]$Fstp[ut] - b[ut])), 0.005)
})
