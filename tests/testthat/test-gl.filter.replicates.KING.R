# gl.filter.replicates.KING removes the samples gl.report.replicates.KING
# identifies, and re-applies other thresholds to the report's pair table.

quiet <- function(expr) {
  capture.output(res <- suppressWarnings(expr))
  res
}

# platypus with two duplicates: T27 copied with 60% heterozygote dropout and
# SUS20 copied exactly
make_dups <- function() {
  x <- platypus.gl
  G <- as.matrix(x)
  set.seed(1)
  g <- G["T27", ]
  h <- which(g == 1)
  k <- sample(h, round(0.6 * length(h)))
  g[k] <- sample(c(0, 2), length(k), replace = TRUE)
  dm <- rbind(G, g, G["SUS20", ])
  rownames(dm) <- c(rownames(G), "T27_dup", "SUS20_dup")
  y <- new("genlight", gen = dm, ploidy = 2)
  indNames(y) <- rownames(dm)
  locNames(y) <- colnames(G)
  pop(y) <- factor(c(as.character(pop(x)), as.character(pop(x)[match(c("T27", "SUS20"), indNames(x))])))
  suppressWarnings(gl.compliance.check(y, verbose = 0))
}

test_that("with the report's thresholds the filter drops exactly ind.list.drop", {
  x <- make_dups()
  rr <- quiet(gl.report.replicates.KING(x, plot.out = FALSE, verbose = 0))
  expect_true(length(rr$ind.list.drop) >= 2)
  y <- quiet(gl.filter.replicates.KING(x, rr, verbose = 0))
  drop <- rr$ind.list.drop
  expect_setequal(setdiff(indNames(x), indNames(y)), drop)
  expect_equal(nInd(y), nInd(x) - length(drop))
  # a single history entry for the call
  expect_equal(length(y@other$history), length(x@other$history) + 1)
  expect_match(deparse(y@other$history[[length(y@other$history)]])[1],
               "gl.filter.replicates.KING")
})

test_that("the filter can loosen and tighten the report's thresholds", {
  x <- make_dups()
  rr <- quiet(gl.report.replicates.KING(x, plot.out = FALSE, verbose = 0))
  tab <- rr$table.all
  # without the heterozygote share only the exact copy is a duplicate
  y <- quiet(gl.filter.replicates.KING(x, rr, het_threshold = 1.01,
                                       verbose = 0))
  expect_true(all(c("T27", "T27_dup") %in% indNames(y)))
  expect_equal(sum(!c("SUS20", "SUS20_dup") %in% indNames(y)), 1)
  # a lower KING threshold flags more pairs, read from table.all
  loose <- 0.1
  fl <- !is.na(tab$KING) & (tab$KING > loose | tab$het_share > 0.8)
  y2 <- quiet(gl.filter.replicates.KING(x, rr, king_threshold = loose,
                                        verbose = 0))
  expect_true(all(setdiff(indNames(x), indNames(y2)) %in%
                    c(tab$ind1[fl], tab$ind2[fl])))
  expect_gte(nInd(x) - nInd(y2), length(rr$ind.list.drop))
})

test_that("the filter returns the object unchanged when nothing is flagged", {
  x <- make_dups()
  rr <- quiet(gl.report.replicates.KING(x, plot.out = FALSE, verbose = 0))
  y <- quiet(gl.filter.replicates.KING(x, rr, king_threshold = 2,
                                       het_threshold = 2, verbose = 0))
  expect_identical(indNames(y), indNames(x))
})

test_that("the filter rejects a report from gl.report.replicates", {
  x <- make_dups()
  old <- quiet(gl.report.replicates(x, plot.out = FALSE, verbose = 0))
  expect_error(quiet(gl.filter.replicates.KING(x, old, verbose = 0)),
               "gl.report.replicates.KING")
})
