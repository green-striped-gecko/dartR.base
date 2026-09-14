#' @name gl.sim.crosses
#' @title Generates crosses between fathers and mothers
#' @family simulation

#' @description
#' Crosses a cohort of fathers (in one genlight object) with a cohort of
#' mothers (in a second genlight object), generating a brood of a specified
#' size for each pair, then retains a random subset of the offspring.

#' @details
#' Pairing is positional: father \emph{i} is crossed with mother \emph{i}, in
#' every brood replicate. The two parent cohorts must therefore hold the same
#' number of individuals, scored on the same loci. Randomising which male
#' breeds with which female is the caller's responsibility, and is done by
#' subsampling the parent cohorts before calling this function.

#' Each parent contributes one gamete per offspring. A homozygous locus
#' transmits its only allele; a heterozygous locus transmits one of its two
#' alleles, drawn independently for every heterozygous call, in every parent,
#' in every brood replicate. A locus scored as missing in either parent yields
#' a missing call in the offspring.

#' This script is to be used in conjunction with gl.subsample.ind() applied
#' initially to a base genlight object containing male and female genotypes.
#' The workflow is:

#' (a) Select the males from the base genlight object using gl.keep.pop() with
#' pop.list set to the label used for males in the sex metadata and the as.pop
#' parameter set to "sex". Note that gl.keep.pop() restores the original
#' population assignments afterwards, and that gl.subsample.ind() subsamples
#' within populations, so assign a single population to the result before
#' subsampling.

#' (b) Select the females in the same way.

#' (c) Subsample a cohort of males and an equally sized cohort of females for
#' breeding using gl.subsample.ind() and the replace parameter as follows.
#' To enforce monogamy -- subsample both cohorts with replace=FALSE.
#' To admit polygyny -- subsample the fathers with replace=TRUE and the
#' mothers with replace=FALSE, so that a male can appear more than once.
#' To admit polyandry -- subsample the fathers with replace=FALSE and the
#' mothers with replace=TRUE.
#' To admit promiscuity -- subsample both cohorts with replace=TRUE.
#' These are simple scenarios that leave the number of mates per individual to
#' chance, depending on the random selection of parents with replacement.

#' (d) Cross the two cohorts using gl.sim.crosses(), retaining a subset of the
#' offspring at random.

#' \preformatted{
#'   males <- gl.keep.pop(testset.gl, pop.list = "Male",
#'                        as.pop = "sex", verbose = 0)
#'   females <- gl.keep.pop(testset.gl, pop.list = "Female",
#'                          as.pop = "sex", verbose = 0)
#'   pop(males) <- rep("males", nInd(males))
#'   pop(females) <- rep("females", nInd(females))
#'   fathers <- gl.subsample.ind(males, n = 10, replace = TRUE, verbose = 0)
#'   mothers <- gl.subsample.ind(females, n = 10, replace = FALSE, verbose = 0)
#'   offspring <- gl.sim.crosses(fathers, mothers, broodsize = 5, n = 20)
#' }

#' The offspring object records its parentage in
#' \code{@@other$ind.metrics$mother} and \code{@@other$ind.metrics$father}.

#' Set error.check to FALSE if using this script in simulations. The structural
#' checks (parameter ranges, matching cohort sizes and locus counts) always
#' run; error.check governs only the comparison of the two locus panels for
#' identity, which scales with the number of loci.

#' @param fathers Genlight object of potential fathers [required].
#' @param mothers Genlight object of potential mothers [required].
#' @param broodsize Number of offspring per mother. Must be a positive integer;
#' any other value is replaced by 10 with a warning [default 10].
#' @param sexratio Expected proportion of female offspring. Must be in the
#' range 0 to 1; any other value is replaced by 0.5 with a warning
#' [default 0.5].
#' @param n Number of offspring to retain, drawn at random from the full brood.
#' If n exceeds the number of offspring generated, all of them are retained,
#' with a warning [default NULL, resolving to 1000 or mothers*broodsize,
#' whichever is the lesser].
#' @param error.check If TRUE, also checks that the two parent cohorts are
#' scored on identical loci [default TRUE].
#' @param compliance.check If TRUE, will perform a compliance check on the
#' resultant genlight object before returning it [default TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @return A genlight object with n offspring of both sexes.

#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @seealso \code{\link{gl.subsample.ind}}, \code{\link{gl.keep.pop}}

#' @examples
#' males <- gl.keep.pop(testset.gl, pop.list = "Male",
#'                      as.pop = "sex", verbose = 0)
#' females <- gl.keep.pop(testset.gl, pop.list = "Female",
#'                        as.pop = "sex", verbose = 0)
#' pop(males) <- rep("males", nInd(males))
#' pop(females) <- rep("females", nInd(females))
#' fathers <- gl.subsample.ind(males, n = 5, replace = FALSE, verbose = 0)
#' mothers <- gl.subsample.ind(females, n = 5, replace = FALSE, verbose = 0)
#' offspring <- gl.sim.crosses(fathers, mothers, broodsize = 4, n = 10,
#'                             verbose = 0)
#' nInd(offspring)
#' head(offspring@other$ind.metrics)

#' @importFrom stats runif
#' @export

gl.sim.crosses <- function(fathers,
                           mothers,
                           broodsize = 10,
                           sexratio = 0.5,
                           n = NULL,
                           error.check = TRUE,
                           compliance.check = TRUE,
                           verbose = NULL) {
  # Preliminaries -------------
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)

  # CHECK DATATYPE
  # Both calls are made for their side effect: they stop on non-SNP input. The
  # algorithm reads a score of 1 as a heterozygote, which has no meaning for
  # presence/absence (SilicoDArT) data scored 0/1 at ploidy 1.
  if (verbose >= 2) {
    cat(report("  father --"))
  }
  utils.check.datatype(fathers, accept = "SNP", verbose = verbose)
  if (verbose >= 2) {
    cat(report("  mother --"))
  }
  utils.check.datatype(mothers, accept = "SNP", verbose = verbose)

  # Function-specific error checking -----------
  # CHECK BROODSIZE
  # The fallback of 10 is assigned, not merely announced: an unvalidated
  # broodsize makes 1:broodsize count down and leaves the offspring unnamed.
  if (length(broodsize) != 1 || is.na(broodsize) ||
      !is.numeric(broodsize) || broodsize < 1 ||
      broodsize != round(broodsize)) {
    if (verbose >= 1) {
      cat(warn("  Warning: Brood size must be a positive integer.",
               "Set to 10\n"))
    }
    broodsize <- 10
  }

  # CHECK SEXRATIO
  if (length(sexratio) != 1 || is.na(sexratio) ||
      !is.numeric(sexratio) || sexratio < 0 || sexratio > 1) {
    if (verbose >= 1) {
      cat(warn("  Warning: Sex ratio must be in the range 0 to 1.",
               "Set to 0.5\n"))
    }
    sexratio <- 0.5
  }

  # CHECK N
  if (!is.null(n)) {
    if (length(n) != 1 || is.na(n) || !is.numeric(n) || n < 1 ||
        n != round(n)) {
      stop(error(
        "Fatal Error: Number of offspring to retain (n) must be a positive",
        "integer\n"
      ))
    }
  }

  # CHECK THAT THE TWO PARENT COHORTS CAN BE PAIRED
  # Crossing is element-wise on the two gamete arrays, so the cohorts must
  # match in size and in the loci they are scored on.
  if (nInd(fathers) != nInd(mothers)) {
    stop(
      error(
        "Fatal Error: Pairing is positional, so the two parent cohorts must",
        "hold the same number of individuals. Found",
        nInd(fathers),
        "fathers and",
        nInd(mothers),
        "mothers\n"
      )
    )
  }
  if (nLoc(fathers) != nLoc(mothers)) {
    stop(
      error(
        "Fatal Error: The two parent cohorts must be scored on the same",
        "loci. Found",
        nLoc(fathers),
        "loci for the fathers and",
        nLoc(mothers),
        "for the mothers\n"
      )
    )
  }
  if (error.check) {
    if (!identical(locNames(fathers), locNames(mothers))) {
      stop(
        error(
          "Fatal Error: The two parent cohorts have the same number of loci",
          "but not the same loci. Check locNames() on both objects\n"
        )
      )
    }
  }

  # CHECK TOTAL NUMBER OF OFFSPRING
  # An unspecified n resolves to the documented lesser of 1000 and the brood
  # total, so it can never exceed what the cross produces.
  noff <- nInd(mothers) * broodsize
  if (is.null(n)) {
    n <- min(1000, noff)
  } else if (noff < n) {
    if (verbose >= 1) {
      cat(warn(
        "  Warning: Sum of broods less than the specified number of",
        "offspring to return, returning",
        noff,
        "\n"
      ))
    }
  }

  # DO THE JOB -------------

  # Draw haploid gametes from a cohort of parents, one gamete per parent per
  # brood replicate. Every heterozygous call gets its own random draw,
  # addressed by position, so no two loci and no two parents share a draw. The
  # array is preallocated rather than grown by rbind() in the loop.
  draw.gametes <- function(pmat, nrep) {
    np <- nrow(pmat)
    gam <- matrix(NA_real_, nrow = np * nrep, ncol = ncol(pmat))
    for (i in seq_len(nrep)) {
      g <- pmat
      het <- which(g == 1)  # which() drops NA, so missing calls stay missing
      g[het] <- sample(c(0, 2), length(het), replace = TRUE)
      gam[((i - 1) * np + 1):(i * np), ] <- g
    }
    gam
  }

  # Generate maternal haplotypes (ova) and paternal haplotypes (sperm)
  ova <- draw.gametes(as.matrix(mothers), broodsize)
  sperm <- draw.gametes(as.matrix(fathers), broodsize)

  # Generate offspring (zygote) genotypes, taking advantage of the dartR coding
  offmat <- (ova + sperm) / 2

  # Record parentage. Gametes are stacked in blocks of one per parent per brood
  # replicate, so row j of block i belongs to parent j.
  mother.id <- rep(indNames(mothers), times = broodsize)
  father.id <- rep(indNames(fathers), times = broodsize)

  # Retain n offspring, drawn at random from the full brood
  if (n < nrow(offmat)) {
    keep <- sort(sample.int(nrow(offmat), n))
    offmat <- offmat[keep, , drop = FALSE]
    mother.id <- mother.id[keep]
    father.id <- father.id[keep]
  }

  # Build the offspring object. Individual names are derived from the matrix
  # itself, so the count cannot disagree with the number of rows.
  ind.names <- paste0("Po_", seq_len(nrow(offmat)))
  gl2 <-
    new(
      "genlight",
      gen = offmat,
      ind.names = ind.names,
      loc.names = locNames(mothers),
      ploidy = rep(2, nrow(offmat))
    )
  if (!is(gl2, "dartR")) {
    class(gl2) <- "dartR"
  }

  # Assign sex. Both levels are declared so that a cohort that happens to be
  # single-sexed still carries a two-level factor.
  sr <- factor(ifelse(runif(nInd(gl2)) < sexratio, "female", "male"),
               levels = c("female", "male"))

  # Populate the metadata. Assigning into @other$ind.metrics before it exists
  # writes into a NULL and yields a bare list, so build the data frame whole.
  gl2@other$ind.metrics <- data.frame(
    id = ind.names,
    sex = sr,
    mother = mother.id,
    father = father.id,
    stringsAsFactors = FALSE
  )
  pop(gl2) <- factor(rep("pop1", nInd(gl2)))

  # Carry the parental locus metrics over, so that sequence-level columns
  # (TrimmedSequence, AlleleID, SnpPosition) survive the rebuild. The
  # frequency-derived columns are stale for the offspring, so all flags are
  # set FALSE and the metrics are recalculated on demand.
  loc.metrics <- mothers@other$loc.metrics
  if (is.null(loc.metrics) || nrow(loc.metrics) != nLoc(gl2)) {
    loc.metrics <- data.frame(row.names = seq_len(nLoc(gl2)))
  }
  gl2@other$loc.metrics <- loc.metrics
  gl2 <- utils.reset.flags(gl2, set = FALSE, verbose = 0)

  if (compliance.check) {
    gl2 <- gl.compliance.check(gl2, verbose = 0)
  }

  # ADD TO HISTORY
  # The object is newly constructed, so its provenance is this call alone. Any
  # entries left by the internal gl.compliance.check() and gl.recalc.metrics()
  # calls reference local variables and are discarded.
  gl2@other$history <- list(match.call())

  # FLAG SCRIPT END ---------------

  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(gl2)
}
