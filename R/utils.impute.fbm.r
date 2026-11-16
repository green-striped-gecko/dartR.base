## Imputation of FBM-backed dartR objects in blocks of loci (used by
## gl.impute). The decoded genotype matrix is never held in memory: each
## pass reads blocks of loci (the storage order of the backing file) and
## writes the imputed values into a copy of the FBM, never into the caller's
## backing file.
##
## Genotype coding: 0, 1, 2 (SNP; number of copies of the alternative
## allele) or 0, 1 (presence-absence), NA missing; FBM code mapping CODE_012,
## where the raw code 3 decodes to NA. Missing values are therefore written
## as 3 (assigning NA to an FBM.code256 would store 0).
##
## Every result is identical to the gen-backed path of gl.impute, including
## the random draws for the same seed:
##   - allele frequencies are integer allele sums divided by the summed
##     ploidy of the scored individuals, as in adegenet::glMean (exact);
##   - random draws are made in the same order as the gen-backed path:
##     column-major within a population for HW (populations in the outer
##     loop, hence one pass per population), column-major over the whole
##     matrix for 'random' and fill.residual. Splitting sample() or rbinom()
##     into consecutive pieces consumes the random stream exactly like one
##     call, so blocking by loci does not change the draws;
##   - nearest-neighbour distances are sums of small integers, which are
##     exact in double precision whatever the block order, so ties (and the
##     neighbour order) are the same as with utils.dist.ind.snp /
##     utils.dist.binary.

## Write a decoded block back into G, missing values as code 3
.impute_write <- function(G, rows, cols, B) {
  B[is.na(B)] <- 3
  if (is.null(rows)) G[, cols] <- B else G[rows, cols] <- B
  invisible(NULL)
}

## Alternative-allele frequency per locus of a block (rows x loci), as
## adegenet::glMean (alleleAsUnit = TRUE): sum of the genotypes divided by
## the summed ploidy of the scored individuals; NaN where none is scored.
.impute_q <- function(B, na, pl_rows) {
  colSums(B, na.rm = TRUE) / colSums((!na) * pl_rows)
}

## Locus (column) of each missing value, in the column-major order of B[na];
## avoids col(B), an integer matrix as large as the block
.impute_na_col <- function(na) {
  (which(na) - 1L) %/% nrow(na) + 1L
}

#' @name utils.impute.fbm.na.by.pop
#' @title Missing data per population of an FBM-backed object, in blocks
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' Counts, in one pass over blocks of loci, the quantities that
#' \code{gl.impute} reports per population, with the same definitions as on
#' a gen-backed object (\code{glNA}, which counts each missing genotype as
#' its individual's ploidy).
#' @param G An FBM.code256 [required].
#' @param groups List of row indices, one element per population
#' [required].
#' @param ploidy Ploidy of each individual (row) [required].
#' @param pl Ploidy assumed by the report: 2 for SNP data, 1 for
#' presence-absence data [required].
#' @return A list with, per population, \code{all_na} (loci where
#' \code{glNA} equals \code{pl} times the number of individuals) and
#' \code{nas} (\code{sum(glNA) / pl}), and \code{total} (the number of
#' missing genotypes in the whole matrix).
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @noRd
utils.impute.fbm.na.by.pop <- function(G, groups, ploidy, pl) {
  n <- nrow(G)
  ng <- length(groups)
  W <- matrix(0, n, ng)                       # ploidy of each row, by group
  for (k in seq_len(ng)) W[groups[[k]], k] <- ploidy[groups[[k]]]
  size_g <- lengths(groups)
  # integer counts as on the gen-backed path (sums of logicals), so that
  # they print the same way
  all_na <- integer(ng)
  nas <- numeric(ng)
  total <- 0L
  for (cols in .fbm_col_blocks(ncol(G), n)) {
    na <- is.na(G[, cols, drop = FALSE])
    total <- total + sum(na)
    M <- crossprod(W, na)                     # groups x loci: glNA per pop
    all_na <- all_na + as.integer(rowSums(M == pl * size_g))
    nas <- nas + rowSums(M)
  }
  list(all_na = all_na, nas = nas / pl, total = total)
}

## The per-population messages that gl.impute prints at verbose >= 2 on the
## gen-backed path, from the counts of utils.impute.fbm.na.by.pop (same text,
## same order).
.impute_fbm_messages <- function(method, pop_names, size_g, st, verbose) {
  if (method == "neighbour") {
    cat(report("  Imputation based on drawing from the nearest neighbour\n"))
  }
  for (k in seq_along(pop_names)) {
    loci_all_nas <- st$all_na[k]
    number_imputations <- st$nas[k] - (loci_all_nas * size_g[k])
    if (method == "frequency") {
      cat(report("  Imputation based on average allele frequencies, population-wise\n"))
    }
    if (method == "HW") {
      cat(report("  Imputation based on average allele HW sampling, population-wise\n"))
    }
    if (loci_all_nas < 1) next
    if (method == "neighbour") {
      cat(warn(
        "  Warning: Population ", pop_names[k],
        " has ", loci_all_nas, " loci with all missing values.\n",
        sep = ""
      ))
      if (verbose >= 3) {
        cat(report(
          "  Method = 'neighbour': ", number_imputations,
          " values to be imputed (excluding all-NA loci).\n",
          sep = ""
        ))
      }
    } else {
      cat(warn("  Warning: Population ", pop_names[k], " has ",
               loci_all_nas, " loci with all missing values.\n"))
      if (verbose >= 3) {
        cat(report(paste0("  Method= '", method, "':"),
                   number_imputations, "values to be imputed.\n"))
      }
    }
  }
  invisible(NULL)
}

#' @name utils.impute.fbm
#' @title Imputes an FBM-backed genotype matrix in blocks of loci
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' The FBM path of \code{\link{gl.impute}} for methods "frequency", "HW",
#' "neighbour" and "random". The imputed genotypes are written into a copy
#' of \code{G}; the result is identical to the gen-backed path, including
#' the random draws for the same seed.
#'
#' Memory: blocks of loci sized by \code{options(bigstatsr.block.sizeGB)}.
#' Method "neighbour" also holds n x n matrices (distances, pairwise counts;
#' 8 bytes per pair each: 32 Mb for 2,000 individuals, 3.2 Gb for 20,000).
#' @param G An FBM.code256 with code mapping CODE_012 [required].
#' @param method One of "frequency", "HW", "neighbour", "random"
#' [required].
#' @param groups List of row indices, one element per population, in the
#' order of \code{seppop} [required].
#' @param ploidy Ploidy of each individual (row) [required].
#' @param datatype "SNP" or "SilicoDArT" [required].
#' @param freq_fill,bernoulli_fill The fill functions of \code{gl.impute}
#' [required].
#' @return A list with \code{G}, the imputed FBM (a new backing file), and
#' \code{remaining}, the number of loci left missing per individual (method
#' "neighbour" only, else NULL).
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @noRd
utils.impute.fbm <- function(G, method, groups, ploidy, datatype,
                             freq_fill, bernoulli_fill) {
  n <- nrow(G)
  p <- ncol(G)
  pl <- if (datatype == "SilicoDArT") 1 else 2
  G2 <- utils.fbm.copy(G)
  remaining <- NULL

  if (method == "frequency") {
    # no random draws: one pass, populations within each block
    for (cols in .fbm_col_blocks(p, n)) {
      B <- G[, cols, drop = FALSE]
      if (!anyNA(B)) next
      for (r in groups) {
        if (!length(r)) next
        Bg <- B[r, , drop = FALSE]
        na <- is.na(Bg)
        if (!any(na)) next
        q <- .impute_q(Bg, na, ploidy[r])
        Bg[na] <- freq_fill(q[.impute_na_col(na)])
        B[r, ] <- Bg
      }
      .impute_write(G2, NULL, cols, B)
    }
  }

  if (method == "HW") {
    # populations in the outer loop: the draws follow the gen-backed order
    # (each population column-major), so one pass per population
    for (r in groups) {
      if (!length(r)) next
      for (cols in .fbm_col_blocks(p, length(r))) {
        Bg <- G[r, cols, drop = FALSE]
        na <- is.na(Bg)
        if (!any(na)) next
        q <- .impute_q(Bg, na, ploidy[r])
        qna <- q[.impute_na_col(na)]
        if (datatype == "SNP") {
          Bg[na] <- unname(unlist(lapply(qna, function(x) {
            sample_genotype(q_freq = x)
          })))
        } else {
          Bg[na] <- bernoulli_fill(qna)
        }
        .impute_write(G2, r, cols, Bg)
      }
    }
  }

  if (method == "random") {
    for (cols in .fbm_col_blocks(p, n)) {
      B <- G[, cols, drop = FALSE]
      na <- is.na(B)
      if (!any(na)) next
      B[na] <- sample(c(0:pl), size = sum(na), replace = TRUE)
      .impute_write(G2, NULL, cols, B)
    }
  }

  if (method == "neighbour") {
    # Euclidean distance over the loci scored in both individuals, as
    # utils.dist.ind.snp and utils.dist.binary (unscaled): with Z the block
    # with NA set to 0 and M the scored indicator,
    #   sumsq = Z^2 M' + M (Z^2)' - 2 Z Z',  L = M M'.
    # Integer sums, exact in doubles, so identical to the gen-backed values.
    S <- matrix(0, n, n)
    L <- matrix(0, n, n)
    for (cols in .fbm_col_blocks(p, n)) {
      Z <- G[, cols, drop = FALSE]
      M <- !is.na(Z)
      Z[!M] <- 0
      storage.mode(M) <- "double"
      Z2 <- Z * Z
      A <- tcrossprod(Z2, M)
      S <- S + A + t(A) - 2 * tcrossprod(Z)
      L <- L + tcrossprod(M)
      rm(Z, M, Z2, A)
    }
    D <- sqrt(S)
    D[L == 0] <- NA
    rm(S, L)
    D[is.na(D)] <- Inf     # undefined distances last, as in gl.impute
    diag(D) <- Inf         # never pick self as neighbour

    # Each locus is filled independently, from the current values of the
    # neighbours at that locus; walking the individuals in order within each
    # block therefore gives the same result as the gen-backed whole-matrix
    # walk (where earlier individuals are already filled).
    remaining <- integer(n)
    ord <- vector("list", n)
    for (cols in .fbm_col_blocks(p, n)) {
      B <- G[, cols, drop = FALSE]
      if (!anyNA(B)) next
      for (i in seq_len(n)) {
        miss <- is.na(B[i, ])
        if (!any(miss)) next
        if (is.null(ord[[i]])) ord[[i]] <- order(D[i, ], decreasing = FALSE)
        for (j in ord[[i]]) {
          if (!any(miss)) break
          cand <- B[j, miss]
          ok <- !is.na(cand)
          if (any(ok)) {
            fill_pos <- which(miss)[ok]
            B[i, fill_pos] <- cand[ok]
            miss[fill_pos] <- FALSE
          }
        }
        remaining[i] <- remaining[i] + sum(miss)
      }
      .impute_write(G2, NULL, cols, B)
    }
  }

  list(G = G2, remaining = remaining)
}

#' @name utils.impute.fbm.residual
#' @title Fills residual missing values of an FBM-backed object
#' @family utilities
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' The \code{fill.residual} step of \code{\link{gl.impute}} in one pass over
#' blocks of loci: each missing value is drawn from the global allele
#' frequency at its locus (\code{s_alleles} for SNP data, a Bernoulli draw
#' with the band frequency for presence-absence data). Loci with no scored
#' genotype stay missing. \code{G} is modified in place.
#' @param G An FBM.code256 with code mapping CODE_012 [required].
#' @param ploidy Ploidy of each individual (row) [required].
#' @param datatype "SNP" or "SilicoDArT" [required].
#' @param bernoulli_fill The fill function of \code{gl.impute} [required].
#' @return \code{G}, invisibly.
#' @author Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @noRd
utils.impute.fbm.residual <- function(G, ploidy, datatype, bernoulli_fill) {
  n <- nrow(G)
  for (cols in .fbm_col_blocks(ncol(G), n)) {
    B <- G[, cols, drop = FALSE]
    na <- is.na(B)
    if (!any(na)) next
    q <- .impute_q(B, na, ploidy)
    qna <- q[.impute_na_col(na)]
    if (datatype == "SNP") {
      B[na] <- unname(unlist(lapply(qna, function(x) {
        as.numeric(s_alleles(q_freq = x))
      })))
    } else {
      B[na] <- bernoulli_fill(qna)
    }
    .impute_write(G, NULL, cols, B)
  }
  invisible(G)
}
