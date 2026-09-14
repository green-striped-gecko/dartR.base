# Characterization tests for gl.tree.fitch (function-review campaign).
# Originally pinned the pre-review behaviour, defects included; updated
# 2026-09-06 when the approved findings of the review were applied — every
# flipped expectation is tagged with the finding it maps to. See
# function-review/reports/dartR.base/gl.tree.fitch.md.
#
# The function drives the external PHYLIP fitch/consense binaries. Tests that
# need the binaries are gated on DARTR_PHYLIP_PATH (fallback: the path used in
# the function's own examples) and skip when absent.
#
# The pre-review file drove the function through a shell() shim because
# system("fitch.exe < fitch.cmd") performed no redirection on Windows; the
# function now uses system2(stdin = ...) and runs unshimmed on all
# platforms. # [approved F1]

phylip_path <- Sys.getenv("DARTR_PHYLIP_PATH", "D:/workspace/R/phylip-3.695/exe")
fitch_bin <- ifelse(.Platform$OS.type == "windows", "fitch.exe", "fitch")
has_phylip <- file.exists(file.path(phylip_path, fitch_bin))

# Deterministic 5-population distance fixture (no missing data), built once
# per file run (gl.dist.phylo dominates the runtime)
.fx_cache <- new.env()
make_fixture <- function() {
  if (is.null(.fx_cache$fx)) {
    gl <- dartR.data::testset.gl
    pops5 <- names(sort(table(adegenet::pop(gl)), decreasing = TRUE))[1:5]
    gl5 <- gl.keep.pop(gl, pop.list = pops5, verbose = 0)
    gl5 <- gl.filter.callrate(gl5, threshold = 1, verbose = 0)
    .fx_cache$fx <- list(gl5 = gl5,
                         D = gl.dist.phylo(gl5, subst.model = "F81",
                                           verbose = 0))
  }
  .fx_cache$fx
}

# Topology pin from the review baseline (PHYLIP fitch 3.695, F81 distances,
# 5 largest testset.gl populations, callrate threshold 1). Branch lengths in
# this string are the pre-review 6-decimal values; only the topology is
# compared (branch lengths are now recovered at higher precision through the
# x10,000 scaling of the distance matrix). # [approved F10]
PINNED_NEWICK <- "(EmmacMDBFo:1e-05,(EmmacMaclG:0,(EmmacBurdM:1e-05,EmmacBurnB:2e-05):5e-05):0,EmmacBrisW:1e-05);"

test_that("non-dist input is rejected by the datatype check", {
  skip_if_not_installed("dartR.data")
  expect_error(
    gl.tree.fitch(D = dartR.data::testset.gl, phylip.path = tempfile(),
                  verbose = 0),
    "expecting dist"
  )
})

test_that("missing PHYLIP binary raises an informative error", {
  skip_if_not_installed("dartR.data")
  # Pre-review: cat(error(...)) then bare stop() left conditionMessage
  # empty; now a proper error condition. # [approved F12]
  d4 <- as.dist(matrix(c(0, 1, 2, 3,
                         1, 0, 4, 5,
                         2, 4, 0, 6,
                         3, 5, 6, 0), 4, 4,
                       dimnames = list(letters[1:4], letters[1:4])))
  err <- tryCatch(
    gl.tree.fitch(D = d4, phylip.path = tempfile(), verbose = 0),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "cannot find")
})

test_that("fewer than four taxa stops with an informative error", {
  # Pre-review: fitch wrote an empty outtree for 3 taxa and the unconditional
  # plot() threw "need finite 'xlim' values"; now guarded before any PHYLIP
  # or plotting work. # [approved F9]
  d3 <- as.dist(matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), 3, 3,
                       dimnames = list(letters[1:3], letters[1:3])))
  expect_error(
    gl.tree.fitch(D = d3, phylip.path = tempfile(), verbose = 0),
    "at least 4 taxa"
  )
  d2 <- as.dist(matrix(c(0, 1, 1, 0), 2, 2,
                       dimnames = list(letters[1:2], letters[1:2])))
  expect_error(
    gl.tree.fitch(D = d2, phylip.path = tempfile(), verbose = 0),
    "at least 4 taxa"
  )
})

test_that("unsupported tree.method and duplicate truncated labels stop", {
  # match.arg rejects anything but FM (pre-review: any string ran FM with an
  # inverted, ungated message). # [approved F14]
  d4 <- as.dist(matrix(1, 4, 4) - diag(4))
  attr(d4, "Labels") <- letters[1:4]
  expect_error(
    gl.tree.fitch(D = d4, phylip.path = tempfile(), tree.method = "NJ",
                  verbose = 0)
  )
  # Labels identical in their first 10 characters are indistinguishable to
  # PHYLIP; pre-review they silently became duplicate tips. # [approved F17]
  d4b <- d4
  attr(d4b, "Labels") <- c("SamePrefix0001", "SamePrefix0002", "c", "d")
  expect_error(
    gl.tree.fitch(D = d4b, phylip.path = tempfile(), verbose = 0),
    "not unique"
  )
})

test_that("end-to-end run: class, tips, pinned topology, out.path, rooting", {
  skip_if_not(has_phylip, "PHYLIP fitch binary not available")
  skip_if_not_installed("dartR.data")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  fx <- make_fixture()

  # Pre-review this call failed on Windows before a tree existed (system()
  # dropped the '< fitch.cmd' redirection); it now completes unshimmed via
  # system2(stdin = ...). # [approved F1]
  tr <- gl.tree.fitch(D = fx$D, phylip.path = phylip_path, verbose = 0)
  expect_s3_class(tr, "phylo")
  expect_setequal(tr$tip.label,
                  trimws(sprintf("%-10s", substr(attr(fx$D, "Labels"), 1, 10))))
  pinned <- ape::read.tree(text = PINNED_NEWICK)
  expect_equal(as.numeric(ape::dist.topo(ape::unroot(tr), ape::unroot(pinned),
                              method = "PH85")), 0)

  # A valid outgroup now roots the returned tree (pre-review pin:
  # is.rooted() was FALSE). # [approved F11]
  tr_og <- gl.tree.fitch(D = fx$D, phylip.path = phylip_path,
                         outgroup = attr(fx$D, "Labels")[2], verbose = 0)
  expect_true(ape::is.rooted(tr_og))

  # out.path now receives the PHYLIP artefacts (pre-review pin: the
  # directory stayed empty). # [approved F6]
  op <- file.path(tempdir(), "gl.tree.fitch-outpath-probe")
  dir.create(op, showWarnings = FALSE)
  invisible(gl.tree.fitch(D = fx$D, phylip.path = phylip_path,
                          out.path = op, verbose = 0))
  expect_true(all(c("infile", "outfile", "outtree") %in% list.files(op)))
})

test_that("bootstrap: supports returned on node.label, cwd restored", {
  skip_if_not(has_phylip, "PHYLIP fitch binary not available")
  skip_if_not_installed("dartR.data")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  fx <- make_fixture()

  wd0 <- getwd()
  set.seed(42)
  tr <- gl.tree.fitch(D = fx$D, x = fx$gl5, bstrap = 3,
                      phylip.path = phylip_path, verbose = 0)

  # The working directory is restored after a bstrap > 1 run (pre-review
  # pin: it was left at tempdir()). # [approved F5]
  expect_identical(normalizePath(getwd()), normalizePath(wd0))

  # Bootstrap supports are attached to the returned tree (pre-review pin:
  # node.label was NULL and the computed values were extracted from the
  # wrong edges). # [approved F2] # [approved F4]
  expect_s3_class(tr, "phylo")
  expect_false(is.null(tr$node.label))
  expect_length(tr$node.label, tr$Nnode)
  sup <- suppressWarnings(as.numeric(tr$node.label))
  expect_true(all(is.na(sup) | (sup >= 0 & sup <= 1)))

  # Independent recomputation of the supports from the replicate trees that
  # fitch wrote (intree, copied to out.path): the proportion of replicate
  # trees containing each bipartition of the returned tree, canonicalised by
  # clade content — must equal node.label exactly. # [approved F2]
  op <- file.path(tempdir(), "gl.tree.fitch-boot-probe")
  dir.create(op, showWarnings = FALSE)
  set.seed(42)
  tr2 <- gl.tree.fitch(D = fx$D, x = fx$gl5, bstrap = 3,
                       phylip.path = phylip_path, out.path = op, verbose = 0)
  expect_identical(tr2$node.label, tr$node.label)
  boot <- ape::read.tree(file.path(op, "intree"))
  splits_of <- function(phy) {
    tips <- phy$tip.label
    lapply(ape::prop.part(phy), function(cl) {
      side <- sort(tips[cl])
      other <- sort(setdiff(tips, side))
      # canonical unrooted split: the lexicographically smaller side
      if (paste(side, collapse = "|") < paste(other, collapse = "|"))
        paste(side, collapse = "|") else paste(other, collapse = "|")
    })
  }
  boot_splits <- lapply(boot, function(b) unique(unlist(splits_of(b))))
  own <- unlist(splits_of(tr2))   # one entry per internal node, in node order
  expected <- vapply(own, function(s)
    sum(vapply(boot_splits, function(bs) s %in% bs, logical(1))),
    numeric(1)) / length(boot)
  got <- as.numeric(tr2$node.label)
  expect_equal(unname(got), unname(expected))

  # The PHYLIP majority-rule consensus tree is attached. # [approved F4]
  expect_s3_class(attr(tr2, "consensus.tree"), "phylo")

  # Returned tree is still the best (non-bootstrap) tree
  pinned <- ape::read.tree(text = PINNED_NEWICK)
  expect_equal(as.numeric(ape::dist.topo(ape::unroot(tr), ape::unroot(pinned),
                              method = "PH85")), 0)
})

test_that("bstrap > 1 without x disables bootstrapping; warning gated", {
  skip_if_not(has_phylip, "PHYLIP fitch binary not available")
  skip_if_not_installed("dartR.data")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  fx <- make_fixture()
  # Pre-review the warning printed even at verbose = 0; it is now gated at
  # verbose >= 1 and verbose = 0 is fully silent. # [approved F7]
  out <- capture.output(
    tr <- gl.tree.fitch(D = fx$D, phylip.path = phylip_path, bstrap = 10,
                        verbose = 0)
  )
  expect_s3_class(tr, "phylo")
  expect_length(out, 0)
  expect_null(tr$node.label)
})

test_that("verbose = 0 is fully silent, including the PHYLIP dialogue", {
  skip_if_not(has_phylip, "PHYLIP fitch binary not available")
  skip_if_not_installed("dartR.data")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  fx <- make_fixture()
  # Pre-review the full PHYLIP menu dialogue printed at every verbosity;
  # subprocess output is now discarded below verbose 2. # [approved F7]
  out <- capture.output(
    tr <- gl.tree.fitch(D = fx$D, x = fx$gl5, bstrap = 3,
                        phylip.path = phylip_path, verbose = 0)
  )
  expect_length(out, 0)
  expect_s3_class(tr, "phylo")
})

test_that("a plotting failure does not destroy the computed tree", {
  skip_if_not(has_phylip, "PHYLIP fitch binary not available")
  skip_if_not_installed("dartR.data")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  fx <- make_fixture()
  # plot.type is invalid, so plot.phylo throws inside the tryCatch; the
  # tree must still be returned (pre-review: an error at the plot step lost
  # the computed tree). # [approved F8]
  tr <- gl.tree.fitch(D = fx$D, phylip.path = phylip_path,
                      plot.type = "no-such-type", verbose = 0)
  expect_s3_class(tr, "phylo")
})

test_that("population mismatch between D and x stops before any PHYLIP run", {
  skip_if_not_installed("dartR.data")
  fx0 <- dartR.data::testset.gl
  d4 <- as.dist(matrix(1, 4, 4) - diag(4))
  attr(d4, "Labels") <- letters[1:4]
  # x's populations are not the taxa of D. # [approved F16]
  expect_error(
    gl.tree.fitch(D = d4, x = fx0, bstrap = 3, phylip.path = tempfile(),
                  verbose = 0),
    "do not match"
  )
})
