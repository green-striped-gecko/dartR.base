#' @name gl.report.shannon
#' @title Reports Shannon-based diversity profiles (Hill numbers) for
#' individuals and populations
#' @family unmatched report

#' @description
#' Calculates diversity of order q (Hill numbers) from SNP dosages, following
#' Ma, Li & Zhang (2020). With level = 'alpha' (the default) the function
#' reports a diversity profile for each individual; with level = 'beta' or
#' 'gamma' it partitions diversity among the individuals within each
#' population and reports one profile per population.

#' @details
#' Genotypes are taken as dosages of the alternate allele (0, 1 or 2), and the
#' dosages act as abundances. For each individual, missing genotypes and zero
#' dosages are dropped, the remaining dosages are normalised to proportions p,
#' and the diversity of order q is computed as
#' \itemize{
#' \item q = 0: the number of loci with non-zero dosage (richness);
#' \item q = 1: exp(-sum(p log p)), the exponential of Shannon entropy
#' (natural logarithms; Hill numbers are free of the logarithm base);
#' \item q >= 2: (sum(p^q))^(1/(1-q)) (q = 2 is the inverse Simpson index).
#' }
#' The profile spans orders q = 0 to order - 1.

#' For level = 'beta' or 'gamma', the individuals of each population form an
#' abundance matrix (individuals x loci, missing genotypes as zero) and the
#' matrix forms of Ma, Li & Zhang (2020) apply per population: gamma is the
#' diversity of the pooled locus abundances, alpha is the mean
#' within-individual diversity (normalised over the population total), and
#' beta = gamma / alpha is the effective number of distinct individuals in the
#' population (ranging from 1, all individuals identical in profile, to the
#' number of individuals). The multiplicative partition gamma = alpha x beta
#' holds at every order. If the object has no population assignments, all
#' individuals are assigned to a single population 'pop1' with a warning.

#' Individuals with no non-missing, non-zero dosages (all genotypes missing or
#' all homozygous for the reference allele) carry no abundance information:
#' their rows are reported as NA at level = 'alpha', they are excluded from the
#' population matrices at level = 'beta' and 'gamma', and a warning names them.

#' If a plot.file is given, the ggplot arising from this function is saved as an
#' "RDS" binary file using saveRDS(); can be reloaded with readRDS(). A file
#' name must be specified for the plot to be saved. If a plot directory
#' (plot.dir) is specified, the ggplot binary is saved to that directory;
#' otherwise to the tempdir().

#' @param x Name of the genlight object containing the SNP data [required].
#' @param plot.display If TRUE, resultant plots are displayed in the plot
#' window [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors Colours for the bars of each diversity order: either a
#' palette function that takes the number of colours, or a vector with one
#' colour per order (order colours, for q = 0 to order - 1)
#' [default gl.colors("dis", verbose = 0)].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL].
#' @param level The level of diversity to report: 'alpha' for a profile per
#' individual; 'beta' or 'gamma' for the partition of diversity among
#' individuals within each population [default 'alpha'].
#' @param order The number of diversity orders to report, spanning q = 0 to
#' order - 1 [default 5].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @return A dataframe, returned invisibly: for level = 'alpha', one row per
#' individual (ID, q0 ... q(order-1)); for level = 'beta' or 'gamma', one row
#' per population (pop, q0 ... q(order-1)).

#' @author Author(s): Ching Ching Lau. Custodian: Ching Ching Lau -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @references
#' \itemize{
#' \item Ma, Z., Li, L., & Zhang, Y. P. (2020). Defining individual-level
#' genetic diversity and similarity profiles. Scientific reports, 10(1), 5805.}

#' @examples
#' require("dartR.data")
#' div <- gl.report.shannon(possums.gl[1:30, ])
#' div.beta <- gl.report.shannon(possums.gl[1:30, ], level = "beta", order = 3)

#' @export

gl.report.shannon <- function(x,
                                plot.display = TRUE,
                                plot.theme = theme_dartR(),
                                plot.colors = gl.colors("dis", verbose = 0),
                                plot.dir = NULL,
                                plot.file = NULL,
                                level = "alpha",
                                order = 5,
                                verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  if (verbose == 0) {
    plot.display <- FALSE
  }

  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir,verbose=0)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  funname <- utils.flag.start(func = funname,
                   build = "Jackson",
                   verbose = verbose)

  # CHECK DATATYPE
  datatype <-
    utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (!is.character(level) || length(level) != 1 ||
      !level %in% c("alpha", "beta", "gamma")) {
    stop(error(
      "  Fatal Error: level must be one of 'alpha', 'beta' or 'gamma'\n"
    ))
  }
  if (!is.numeric(order) || length(order) != 1 || is.na(order) ||
      order < 1 || order != round(order)) {
    stop(error(
      "  Fatal Error: order must be a positive whole number; the profile spans orders q = 0 to order - 1\n"
    ))
  }
  # One colour per diversity order (q = 0 to order - 1)
  if (is.function(plot.colors)) {
    plot.colors <- plot.colors(order)
  }
  if (length(plot.colors) < order) {
    stop(error(
      "  Fatal Error: plot.colors must supply at least", order,
      "colours, one per diversity order (q = 0 to order - 1)\n"
    ))
  }

  # For the population-level partition, population assignments are required
  if (level %in% c("beta", "gamma")) {
    if (is.null(pop(x)) || length(pop(x)) != nInd(x) || any(is.na(pop(x)))) {
      if (verbose >= 2) {
        cat(
          warn(
            "  No population assignments detected, individuals assigned to a single population labelled 'pop1'\n"
          )
        )
      }
      pop(x) <- factor(array("pop1", dim = nInd(x)))
    }
  }

  # Codes of calculation are from Online Supplementary Information (OSI) for:
  # Ma ZS, Li LW and Zhang YP (2019) Defining Individual-Level Genetic
  # Diversity and Similarity Profiles. Scientific Reports

  # Hill numbers (alpha, beta = gamma / alpha, gamma) from abundances given
  # as their distinct values 'v' (> 0) with multiplicities 'w', so that a
  # sum over the abundance vector, sum(f(a)), is sum(w * f(v)). Dosages are
  # 0/1/2 (missing as zero), so an individual's abundances are n1 ones and
  # n2 twos, and a population's pooled abundances (colSums) are integers
  # 1..2N. Same formulas as the OSI code (previously d.chao() on the
  # decoded vector or matrix); only the order of summation differs
  # (relative differences of about 1e-14).
  #   tot: total abundance; N: number of individuals (rows of A);
  #   v, w: element abundances (alpha); vc, wc: pooled abundances (gamma)
  d.chao.counts <- function(tot, N, v, w, vc = NULL, wc = NULL, lev, q) {
    if (lev == 'alpha') {
      keep <- w > 0
      eA <- v[keep] / tot
      w <- w[keep]
      if (q != 1) {
        D.value <- (1 / N) * (sum(w * eA ^ q)) ^ (1 / (1 - q))
      } else{
        D.value <- exp(-sum(w * eA * log(eA)) - log(N))
      }
    }
    if (lev == 'beta') {
      D.value <- d.chao.counts(tot, N, v, w, vc, wc, lev = 'gamma', q) /
        d.chao.counts(tot, N, v, w, vc, wc, lev = 'alpha', q)
    }
    if (lev == 'gamma') {
      keep <- wc > 0 & vc > 0
      ecA <- vc[keep] / tot
      wc <- wc[keep]
      if (q != 1) {
        D.value <- (sum(wc * ecA ^ q)) ^ (1 / (1 - q))
      } else{
        D.value <- exp(-sum(wc * ecA * log(ecA)))
      }
    }
    D.value
  }

  # Per individual, the numbers of dosages 1 and 2, summed over blocks of
  # loci (.gl_rowwise), so a file-backed object is not decoded as a whole
  # (previously the whole genotype matrix was decoded with as.matrix())
  n12 <- .gl_rowwise(x, function(t) {
    cbind(rowSums(t == 1, na.rm = TRUE), rowSums(t == 2, na.rm = TRUE))
  })
  ID <- indNames(x)
  list_order <- order-1

  # Individuals with no non-missing, non-zero dosages carry no abundance
  # information: diversity is undefined for them
  informative <- (n12[, 1] + n12[, 2]) > 0
  if (any(!informative) && verbose >= 1) {
    cat(
      warn(
        "  Warning: individual(s)",
        paste(ID[!informative], collapse = ", "),
        "have no non-missing, non-zero dosages; diversity is undefined and reported as NA\n"
      )
    )
  }

  if (level == "alpha") {
    # SNP diversity per individual
    if (verbose >= 2) {
      cat(
        report(
          "  Calculating SNP diversity for each individual\n"
        )
      )
    }

    div_mat <- matrix(0, nInd(x), order)
    for (n in seq_len(nInd(x))) {
      # the individual's abundances: n1 dosages of 1 and n2 of 2 (N = 1)
      w <- n12[n, ]
      for (q in 0:list_order) {
        div_mat[n, q + 1] <- d.chao.counts(tot = sum(w * c(1, 2)), N = 1,
                                           v = c(1, 2), w = w,
                                           lev = level, q = q)
      }
    }
    div_mat[!informative, ] <- NA

    # output
    div_mat <- as.data.frame(div_mat)
    div_mat <- as.data.frame(cbind(ID, div_mat))
    colnames(div_mat) <- c("ID", paste0("q", c(0:list_order)))

  } else {
    # Partition of diversity among individuals within each population:
    # each population's individuals form an abundance matrix
    # (individuals x loci, missing genotypes as zero) and the matrix forms
    # of Ma, Li & Zhang (2020) give gamma (pooled), alpha (mean
    # within-individual) and beta = gamma / alpha per population
    if (verbose >= 2) {
      cat(
        report(
          "  Calculating", level,
          "diversity across individuals within each population\n"
        )
      )
    }

    pops <- levels(pop(x))
    pop_vec <- as.character(pop(x))
    # Per population, the distribution of the pooled abundances (colSums
    # of the informative individuals, missing genotypes as zero) over loci:
    # tabulated in one pass over blocks of loci and added, so neither the
    # object nor a population is decoded as a whole
    pop_keep <- lapply(pops, function(p) which(pop_vec == p & informative))
    hist_blocks <- .gl_colblocks(x, function(B) {
      lapply(pop_keep, function(r) {
        A <- B[r, , drop = FALSE]
        A[is.na(A)] <- 0
        tabulate(colSums(A) + 1, nbins = 2 * length(r) + 1)
      })
    })
    pop_hist <- lapply(seq_along(pops), function(g) {
      Reduce(`+`, lapply(hist_blocks, function(b) b[[g]]))
    })
    rm(hist_blocks)
    div_mat <- matrix(NA_real_, length(pops), order)
    for (n in seq_along(pops)) {
      keep <- pop_vec == pops[n] & informative
      if (!any(keep)) {
        if (verbose >= 1) {
          cat(
            warn(
              "  Warning: population", pops[n],
              "has no informative individuals; reported as NA\n"
            )
          )
        }
        next
      }
      # element abundances: the dosages 1 and 2 of the population's
      # informative individuals; pooled abundances: values 0..2N with
      # their numbers of loci
      w <- colSums(n12[keep, , drop = FALSE])
      vc <- seq_along(pop_hist[[n]]) - 1
      for (q in 0:list_order) {
        div_mat[n, q + 1] <- d.chao.counts(tot = sum(w * c(1, 2)),
                                           N = sum(keep),
                                           v = c(1, 2), w = w,
                                           vc = vc, wc = pop_hist[[n]],
                                           lev = level, q = q)
      }
    }

    # output
    div_mat <- as.data.frame(div_mat)
    div_mat <- as.data.frame(cbind(pop = pops, div_mat))
    colnames(div_mat) <- c("pop", paste0("q", c(0:list_order)))
  }

  div_mat2 <- reshape2::melt(div_mat, id.vars = colnames(div_mat)[1])
  Ord <- SNP_diversity <- ID_col <- NA
  colnames(div_mat2) <- c("ID_col", "Ord", "SNP_diversity")
  unit <- ifelse(level == "alpha", "individual", "population")

  p1 <-
    ggplot(div_mat2, aes(
      x = ID_col,
    )) + geom_bar(aes(y = SNP_diversity, fill=Ord), position = "dodge",
                  stat = "identity") + plot.theme + theme(
                    axis.ticks.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.ticks.y = element_blank(),
                    legend.position = "none"
                  ) + facet_grid(~Ord) +
    scale_fill_manual(values = plot.colors) +
    labs(fill = "Order", y=paste0(level, " diversity")) +
    ggtitle(paste0(level, " diversity per ", unit, " for different orders"))


  # Optionally save the plot ---------------------

  if(!is.null(plot.file)){
    tmp <- utils.plot.save(p1,
                           dir=plot.dir,
                           file=plot.file,
                           verbose=verbose)
  }

  if (verbose >= 3) {
    cat(report("  Returning a dataframe with SNP diversity values\n"))}

  # PRINTING OUTPUTS
  if (plot.display) {
    suppressWarnings(print(p1))}

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN
  return(invisible(div_mat))}
