# Characterization tests for gl.dist.phylo
# Baseline snapshotted before review (phylogeny chain review, dev at
# ddaed27). Pins current behaviour, defects included; assertions tagged
# [pins defect] are expected to flip when the matching approved finding
# is applied. Anchor values were computed on this machine at the
# reviewed state and are 8 dp rounded (distances are of order 1e-3).
#
# Pipeline pinned here: gl2fasta(method = 1) assembles one concatenated
# sequence per individual from TrimmedSequence + SnpPosition (0-based),
# heterozygotes as IUPAC ambiguity codes, NA genotypes as an all-N tag;
# ape::dist.dna computes the distances (ambiguity codes are treated as
# missing data by ape); by.pop = TRUE averages the individual distances
# between populations, labelled by popNames().

# Deterministic 12-individual, 35-locus platypus fixture: no missing
# data, 3 populations x 4 individuals, SNP position strictly inside the
# tag for every locus (so the overshoot pre-filter removes nothing).
fixA <- function() {
  p <- gl.filter.callrate(platypus.gl, threshold = 1, verbose = 0)
  p <- gl.filter.monomorphs(p, verbose = 0)
  lm <- p@other$loc.metrics
  ok <- which(as.numeric(as.character(lm$SnpPosition)) + 1 <=
                nchar(as.character(lm$TrimmedSequence)))
  keep <- ok[1:40]
  pA <- p[, keep]
  pA@other$loc.metrics <- p@other$loc.metrics[keep, , drop = FALSE]
  sel <- unlist(lapply(split(seq_len(nInd(pA)), pop(pA)), head, 4))
  pA <- gl.keep.ind(pA, ind.list = indNames(pA)[sel], verbose = 0)
  gl.filter.monomorphs(pA, verbose = 0)
}

test_that("population-level anchors per substitution model", {
  pA <- fixA()
  expect_equal(nInd(pA), 12)
  expect_equal(nLoc(pA), 35)
  anchors <- list(
    F81 = c(0.00077072, 0.00107468, 0.00066111),
    K80 = c(0.00077086, 0.00107508, 0.00066132),
    raw = c(0.00077015, 0.00107366, 0.00066066)
  )
  for (m in names(anchors)) {
    D <- gl.dist.phylo(pA, subst.model = m, verbose = 0)
    expect_s3_class(D, "dist", exact = FALSE)
    expect_equal(attr(D, "Size"), 3, info = m)
    expect_identical(labels(D),
                     c("SEVERN_ABOVE", "SEVERN_BELOW", "TENTERFIELD"),
                     info = m)
    expect_equal(round(as.vector(D), 8), anchors[[m]], info = m)
  }
  # deterministic: method 1 assembly has no RNG
  D1 <- gl.dist.phylo(pA, subst.model = "F81", verbose = 0)
  D2 <- gl.dist.phylo(pA, subst.model = "F81", verbose = 0)
  expect_identical(as.vector(D1), as.vector(D2))
})

test_that("individual-level anchors, labels and symmetry", {
  pA <- fixA()
  D <- gl.dist.phylo(pA, subst.model = "F81", by.pop = FALSE, verbose = 0)
  expect_s3_class(D, "dist", exact = FALSE)
  expect_equal(attr(D, "Size"), 12)
  # labels are indName_pop composites, in individual order, not
  # indNames(x) — undocumented contract, pinned as is
  expect_identical(labels(D)[1:3],
                   c("T27_TENTERFIELD", "T35_TENTERFIELD",
                     "SDS4_SEVERN_BELOW"))
  expect_equal(round(as.vector(D)[1:6], 8),
               c(0.00132334, 0.00087925, 0.00132159, 0.00088119,
                 0.00132392, 0.00176263))
  expect_true(isSymmetric(as.matrix(D)))
  expect_true(all(diag(as.matrix(D)) == 0))
  # K80 individual anchors
  DK <- gl.dist.phylo(pA, subst.model = "K80", by.pop = FALSE, verbose = 0)
  expect_equal(round(as.vector(DK)[1:6], 8),
               c(0.00132392, 0.00087951, 0.00132217, 0.00088120,
                 0.00132450, 0.00176367))
})

test_that("population labels track non-alphabetical factor levels", {
  pA <- fixA()
  Dfwd <- gl.dist.phylo(pA, subst.model = "F81", verbose = 0)
  pR <- pA
  pop(pR) <- factor(as.character(pop(pA)),
                    levels = c("TENTERFIELD", "SEVERN_BELOW",
                               "SEVERN_ABOVE"))
  DR <- gl.dist.phylo(pR, subst.model = "F81", verbose = 0)
  expect_identical(labels(DR),
                   c("TENTERFIELD", "SEVERN_BELOW", "SEVERN_ABOVE"))
  mF <- as.matrix(Dfwd)
  mR <- as.matrix(DR)
  # same values under name lookup regardless of level order
  expect_equal(mR[rownames(mF), colnames(mF)], mF)
})

test_that("pairwise.missing anchors on a deterministic-NA fixture", {
  pB <- fixA()
  gmx <- as.matrix(pB)
  idx <- cbind(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3),
               c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 2, 5, 8))
  gmx[idx] <- NA
  pB@gen <- methods::new("genlight", gmx, ploidy = 2)@gen
  Dt <- gl.dist.phylo(pB, subst.model = "F81", pairwise.missing = TRUE,
                      verbose = 0)
  expect_equal(round(as.vector(Dt), 8),
               c(0.00074003, 0.00104375, 0.00066374))
  # [pins defect] complete deletion drops every site at which ANY
  # individual carries an ambiguity code (heterozygote) or N: on this
  # fixture no variable site survives and all distances collapse to 0
  Df <- gl.dist.phylo(pB, subst.model = "F81", pairwise.missing = FALSE,
                      verbose = 0)
  expect_equal(as.vector(Df), c(0, 0, 0))
})

test_that("heterozygous sites contribute nothing to the distance", {
  # [pins defect] method-1 ambiguity codes are treated as missing data
  # by ape::dist.dna, so an individual heterozygous at 10 loci is at
  # distance 0 from one homozygous-reference at all loci
  pA <- fixA()
  gm2 <- as.matrix(pA)[1:2, , drop = FALSE]
  gm2[1, ] <- 0
  gm2[2, ] <- 0
  gm2[2, 1:10] <- 1
  p2 <- pA[1:2, ]
  p2@other$loc.metrics <- pA@other$loc.metrics
  p2@gen <- methods::new("genlight", gm2, ploidy = 2)@gen
  Dhet <- gl.dist.phylo(p2, subst.model = "raw", by.pop = FALSE,
                        verbose = 0)
  expect_equal(as.numeric(Dhet), 0)
  # the same 10 loci as homozygous-alternate register fully
  gm2[2, 1:10] <- 2
  p2@gen <- methods::new("genlight", gm2, ploidy = 2)@gen
  Dhom <- gl.dist.phylo(p2, subst.model = "raw", by.pop = FALSE,
                        verbose = 0)
  expect_equal(round(as.numeric(Dhom), 8), 0.00436872)
})

test_that("min.tag.len filters short tags before assembly", {
  pA <- fixA()
  D <- gl.dist.phylo(pA, subst.model = "F81", min.tag.len = 40,
                     verbose = 0)
  expect_s3_class(D, "dist", exact = FALSE)
  expect_equal(round(as.vector(D), 8),
               c(0.00078193, 0.00109033, 0.00067082))
})

test_that("structure and error paths", {
  pA <- fixA()
  # BH87 at individual level returns ape's asymmetric matrix, not the
  # documented dist object [pins defect]
  Dbh <- gl.dist.phylo(pA, subst.model = "BH87", by.pop = FALSE,
                       verbose = 0)
  expect_true(is.matrix(Dbh))
  # invalid model errors from ape with its model list
  expect_error(
    gl.dist.phylo(pA, subst.model = "ugpma", by.pop = FALSE, verbose = 0),
    "'model' must be one of"
  )
  # [approved F3] SilicoDArT is rejected early by this function's own
  # (repaired) branch with the redirect message, before any processing
  wd0 <- getwd()
  on.exit(setwd(wd0), add = TRUE)
  expect_error(gl.dist.phylo(testset.gs, verbose = 0),
               "works only with SNP data")
  # [approved F2] no setwd remains, so the failure leaves the working
  # directory untouched
  expect_true(identical(getwd(), wd0))
  # [approved F4] a single population with by.pop = TRUE fails with an
  # informative message naming the requirement
  p1 <- gl.keep.pop(pA, pop.list = popNames(pA)[1], verbose = 0)
  expect_error(gl.dist.phylo(p1, verbose = 0),
               "at least two populations")
  expect_true(identical(getwd(), wd0))
  # missing TrimmedSequence column: fails mid-pipeline inside gl2fasta;
  # [approved F2] the working directory survives the mid-pipeline failure
  pnt <- pA
  pnt@other$loc.metrics$TrimmedSequence <- NULL
  expect_error(gl.dist.phylo(pnt, verbose = 0),
               "must include Trimmed Sequences")
  expect_true(identical(getwd(), wd0))
})

test_that("all-monomorphic input fails with an informative error", {
  # [approved F5] the monomorph probe no longer calls
  # gl.filter.monomorphs (which cannot return a zero-locus object); an
  # all-monomorphic object now gets a diagnosis naming the cause
  p <- gl.filter.callrate(platypus.gl, threshold = 1, verbose = 0)
  gm <- as.matrix(p)
  mono <- which(apply(gm, 2, function(cc) {
    length(unique(cc[!is.na(cc)])) == 1
  }))[1:10]
  pm <- p[, mono]
  pm@other$loc.metrics <- p@other$loc.metrics[mono, , drop = FALSE]
  pm <- gl.keep.ind(pm, ind.list = indNames(pm)[1:6], verbose = 0)
  expect_error(gl.dist.phylo(pm, verbose = 0),
               "no polymorphic loci")
})

test_that("verbose = 0 is silent and the return is visible", {
  pA <- fixA()
  # suppressWarnings covers capture.output's own cleanup: gl2fasta's
  # on.exit sink guard removes the caller's sink (defect noted in the
  # review, owned by gl2fasta)
  out <- suppressWarnings(capture.output(
    res <- gl.dist.phylo(pA, verbose = 0)
  ))
  # gl.filter.overshoot at ddaed27 prints "There were no loci with SNP
  # falling outside the trimmed sequence" ungated (its own defect, fixed
  # in the pending overshoot review); filter that known third-party leak
  # -- gl.dist.phylo itself must contribute no output at verbose = 0
  out <- out[!grepl("falling outside the trimmed sequence", out)]
  expect_length(out, 0)
  expect_s3_class(res, "dist", exact = FALSE)
  v <- withVisible(gl.dist.phylo(pA, verbose = 0))
  expect_true(v$visible)
})

test_that("het-load warning and gamma/variance pass-through", {
  pA <- fixA()
  gm <- as.matrix(pA)
  n.het <- sum(gm == 1, na.rm = TRUE)
  n.called <- sum(!is.na(gm))
  pct <- round(100 * n.het / n.called, 1)
  # [approved F1] a warning reports the fraction of heterozygous calls
  # at verbose >= 1 (printed before gl2fasta, so capture.output holds it)
  out1 <- suppressWarnings(capture.output(
    d <- gl.dist.phylo(pA, verbose = 1)
  ))
  expect_true(any(grepl(
    paste0(n.het, " of ", n.called, " genotype calls \\(", pct, "%\\)"),
    out1
  )))
  expect_false(any(grepl("deleted for ALL individuals", out1)))
  # [approved F1] the stronger warning fires when pairwise.missing =
  # FALSE makes the deletion global
  out2 <- suppressWarnings(capture.output(
    d <- gl.dist.phylo(pA, pairwise.missing = FALSE, verbose = 1)
  ))
  expect_true(any(grepl("deleted for ALL individuals", out2)))
  # [approved F6] gamma and variance now exist and reach ape::dist.dna:
  # a gamma-corrected K80 distance differs from the uncorrected one
  Dg <- gl.dist.phylo(pA, subst.model = "K80", gamma = 0.5,
                      by.pop = FALSE, verbose = 0)
  D0 <- gl.dist.phylo(pA, subst.model = "K80",
                      by.pop = FALSE, verbose = 0)
  # zero distances stay zero under the correction; all others inflate
  expect_true(all(as.vector(Dg) >= as.vector(D0)))
  expect_true(any(as.vector(Dg) > as.vector(D0)))
  # variance = TRUE attaches ape's variance attribute at by.pop = FALSE
  Dv <- gl.dist.phylo(pA, subst.model = "K80", variance = TRUE,
                      by.pop = FALSE, verbose = 0)
  expect_false(is.null(attr(Dv, "variance")))
})
