#' @name gl.sim.genotypes
#' @title Generate random genotypes
#' @family data manipulation

#' @description
#' Generates random genotypes for each population of a genlight object by
#' drawing, locus by locus, from the allele frequencies of that population.

#' @details
#' Allele frequencies are estimated WITHIN each population of x, and each
#' population is then simulated separately. For every locus, two haplotype
#' vectors of length n.ind are drawn independently, each with probability
#' equal to the allele frequency of that population at that locus, and the
#' two are added to give a dosage of 0, 1 or 2. The simulated genotypes are
#' therefore in Hardy-Weinberg equilibrium within each population, and the
#' allele frequency differences among the populations of x -- the population
#' structure -- are carried through to the simulated object.

#' n.ind is the number of individuals simulated PER POPULATION, so an object
#' with nPop(x) populations returns n.ind * nPop(x) individuals. The simulated
#' individuals are named Ind_1 ... Ind_n and are assigned the population names
#' of the source object. Locus names, and the number of loci, are those of the
#' source object.

#' A locus with no called genotypes in a given population has no estimable
#' allele frequency. Such a locus is returned as missing (NA) for every
#' individual of that population rather than aborting the run, so a locus that
#' is all-missing across the whole source object is all-missing in the
#' simulated object. No other missing data are simulated: every locus with an
#' estimable frequency is called for every individual of that population.

#' If n.ind exceeds the number of loci, n.ind is reduced to the number of loci.
#' Ordinations such as PCA require more attributes (loci) than entities
#' (individuals), and the cap keeps each simulated population inside that
#' bound.

#' The returned object carries the call that made it as its only history
#' entry, so gl.print.history() on a saved simulated dataset reports the
#' simulation rather than the internal housekeeping calls.

#' @param x Name of the genlight object [required].
#' @param n.ind Number of individuals to be simulated in each population. Must
#' be a single whole number of 1 or more; a non-integer value is rounded with
#' a warning. Capped at the number of loci [default 200].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @return A genlight object holding the simulated genotypes, with n.ind
#' individuals in each population of the source object.

#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @seealso \code{\link{gl.sim.crosses}}, \code{\link{gl.filter.allna}}

#' @examples
#' # possums.gl has 10 populations, so this returns 10 x 20 = 200 individuals
#' sim <- gl.sim.genotypes(possums.gl, n.ind = 20, verbose = 0)
#' nInd(sim)
#' table(pop(sim))
#' # Loci with no calls in a population come back missing for that population
#' sim2 <- gl.sim.genotypes(testset.gl, n.ind = 10, verbose = 0)
#' nLoc(sim2) == nLoc(testset.gl)

#' @export

gl.sim.genotypes <- function(x,
                             n.ind = 200,
                             verbose = NULL) {
  # Preliminaries -------------
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   verbose = verbose)

  # CHECK DATATYPE
  # The algorithm is diploid dosage arithmetic: two haplotypes are added to
  # give a score of 0, 1 or 2 at ploidy 2. That has no meaning for
  # presence/absence (SilicoDArT) data scored 0/1 at ploidy 1.
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # Function-specific error checking -----------
  # CHECK N.IND
  # Unvalidated values reach base R and fail there with messages that name
  # neither n.ind nor this function ("subscript out of bounds", "negative
  # length vectors are not allowed"), or are silently truncated.
  if (is.null(n.ind) || length(n.ind) != 1 || !is.numeric(n.ind) ||
      !is.finite(n.ind)) {
    stop(error(
      "Fatal Error: n.ind must be a single finite numeric value\n"
    ))
  }
  if (n.ind < 1) {
    stop(error(
      "Fatal Error: n.ind must be 1 or greater. Set to", n.ind, "\n"
    ))
  }
  if (n.ind != round(n.ind)) {
    if (verbose >= 1) {
      cat(warn("  Warning: n.ind must be a whole number. Rounding",
               n.ind, "to", round(n.ind), "\n"))
    }
    n.ind <- round(n.ind)
  }

  # Bring the input object up to dartR compliance. A genlight not built by
  # dartR has no loc.metrics.flags, and gl.allele.freq() fails on it with
  # "argument is of length zero".
  x <- gl.compliance.check(x, verbose = 0)

  n.loc <- nLoc(x)

  # CHECK N.IND AGAINST N.LOC
  # The cap is applied, not merely announced.
  if (n.ind > n.loc) {
    if (verbose >= 1) {
      cat(warn("  Warning: the number of individuals (entities) exceeds the number of loci (attributes) which will cause issues for analyses like PCA\n"))
      cat(warn("    Setting n.ind to", n.loc, "\n"))
    }
    n.ind <- n.loc
  }

  # DO THE JOB --------------

  # Extract the allele frequencies for each population separately. Pooling
  # across populations would replace the structure in x with a single
  # panmictic gene pool, and so erase any Wahlund deficit the source carries.
  m <- gl.allele.freq(x, by = "popxloc", verbose = 0)
  m$popn <- as.character(m$popn)
  pop.names <- popNames(x)
  pop.names <- pop.names[pop.names %in% m$popn]
  n.pop <- length(pop.names)

  # One frequency vector per population, in locus order
  freq <- lapply(pop.names, function(p) {
    mp <- m[m$popn == p, , drop = FALSE]
    mp[order(mp$loc_order), "frequency"]
  })
  names(freq) <- pop.names

  # Create an array to hold the new genotypes. It is initialised to NA, so a
  # locus with no estimable frequency in a population is left missing.
  v <- array(NA_real_, dim = c(n.loc, n.ind * n.pop))

  # Populate the array, one population at a time
  for (k in seq_len(n.pop)) {
    f <- freq[[k]]
    cols <- ((k - 1) * n.ind + 1):(k * n.ind)
    for (i in 1:n.loc) {
      # A locus with no calls in this population has no allele frequency to
      # sample from. Leave it missing.
      if (is.na(f[i])) {
        next
      }
      # Generate a vector of individuals, sampling locus i using the allele
      # frequencies for that locus in this population
      v1 <- sample(c(0, 1), size = n.ind, replace = TRUE,
                   prob = c((1 - f[i]), f[i]))
      # Sample again
      v2 <- sample(c(0, 1), size = n.ind, replace = TRUE,
                   prob = c((1 - f[i]), f[i]))
      # Combine the two random haplotypes to form a genotype
      v[i, cols] <- v1 + v2
    }
  }

  # Convert the array v to a new genlight object. Ploidy is per individual.
  ind.names <- paste0("Ind_", 1:(n.ind * n.pop))
  pop.vec <- rep(pop.names, each = n.ind)

  gl <- new(
    "genlight",
    gen = t(v),
    ind.names = ind.names,
    loc.names = locNames(x),
    ploidy = rep(2, n.ind * n.pop)
  )
  if (!is(gl, "dartR")) {
    class(gl) <- "dartR"
  }

  # Populate the metadata. Assigning into @other$ind.metrics before it exists
  # writes into a NULL and yields a bare list, so build the data frame whole.
  gl@other$ind.metrics <- data.frame(
    id = ind.names,
    pop = pop.vec,
    stringsAsFactors = FALSE
  )
  pop(gl) <- factor(pop.vec, levels = pop.names)

  # An empty-but-shaped locus metrics table, one row per locus. Left NULL, it
  # is replaced downstream by a frame whose only column is named
  # "array(NA, nLoc(x))", with a matching "array(NA, 1)" in the flags table.
  # Every metric describes the simulated data, so all flags are set FALSE and
  # the metrics are recalculated by the compliance check.
  gl@other$loc.metrics <- data.frame(row.names = seq_len(n.loc))
  gl@other$loc.metrics.flags <- data.frame(row.names = 1L)
  gl <- utils.reset.flags(gl, set = FALSE, verbose = 0)

  # Enforce its compliance with dartR
  gl <- gl.compliance.check(gl, verbose = 0)

  # ADD TO HISTORY
  # The object is newly constructed, so its provenance is this call alone. The
  # entries left by the internal gl.compliance.check() and gl.recalc.metrics()
  # calls reference local variables and are discarded.
  gl@other$history <- list(match.call())

  # Results summary ------------
  if (verbose >= 3) {
    mat <- as.matrix(gl)
    n.na.pl <- sum(is.na(mat))
    ho <- mean(colMeans(mat == 1, na.rm = TRUE), na.rm = TRUE)
    cat(report("  Simulated", n.ind, "individuals in each of", n.pop,
               "population(s):", nInd(gl), "individuals x", n.loc,
               "loci\n"))
    cat(report("  Allele frequencies estimated from", nInd(x),
               "source individuals\n"))
    cat(report("  Genotypes left missing (locus not called in the source",
               "population):", n.na.pl, "of", length(mat), "\n"))
    cat(report("  Mean simulated heterozygosity:", round(ho, 4), "\n"))
  }

  # FLAG SCRIPT END ---------------

  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(gl)

}
