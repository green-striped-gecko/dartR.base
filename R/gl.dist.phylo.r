#' @name gl.dist.phylo
#' @title Generates a distance matrix from a SNP genlight object taking into
#' account a substitution model
#'
#' @family phylogeny
#'
#' @description
#' Generates a distance matrix for individuals or populations in a genlight
#' object using one of a selection of substitution models.
#'
#' @details
#' The script takes a genlight object as input, creates a set of sequences
#' from the trimmed sequence tags for each individual, calculates distances
#' between the individuals and then optionally averages those distances between
#' the populations defined in the genlight object (typically OTUs).
#'
#' The input genlight object must carry the locus metrics 'TrimmedSequence'
#' and 'SnpPosition' in the @other$loc.metrics slot, and the allele pairs
#' (e.g. 'G/A') in the @loc.all slot. Objects read from DArT reports have
#' these; genlight objects from other sources (for example, from a VCF file)
#' may not, in which case this function cannot be used. Note also that
#' secondaries (multiple SNPs on one sequence tag) each contribute a full
#' duplicate copy of the tag to the assembled sequences, with only their own
#' SNP substituted; consider gl.filter.secondaries() before calling this
#' function.
#'
#' Heterozygous genotypes are written into the assembled sequences as IUPAC
#' ambiguity codes, and ape::dist.dna() treats ambiguity codes as missing
#' data. Heterozygous sites therefore carry no distance signal -- distances
#' are driven by homozygous differences only, and are biased downward
#' wherever heterozygosity is appreciable. With pairwise.missing = FALSE the
#' effect is global: any site at which any individual is heterozygous (or
#' missing) is deleted for all individuals, which can remove most or all of
#' the variable sites. A warning reporting the fraction of heterozygous
#' genotype calls is printed at verbose 1 or above.
#'
#' min.tag.len : Sequence tags can vary considerably in length, which
#' results in large numbers of Ns in alignments. This can have an impact of
#' distance measures depending on how missing values are managed. To minimize
#' this effect, you might elect to filter on tag length using this parameter.
#'
#' subst.model : Use this parameter to specify the substitution model,
#' selecting from the list used by package ape.
#'
#'      \itemize{
#'      \item
#'       raw: This is simply the proportion or the number of sites that differ
#'       between each pair of sequences. This may be useful to draw "saturation
#'       plots". The options variance and gamma have no effect, but
#'       pairwise.missing can.
#'      \item
#'       TS, TV: These are the numbers of transitions and transversions,
#'       respectively.
#'      \item
#'       JC69: This model was developed by Jukes and Cantor (1969). It assumes
#'       that all substitutions (i.e. a change of a base by another one) have
#'       the same probability. This probability is the same for all sites along
#'       the DNA sequence. This last assumption can be relaxed by assuming that
#'       the substition rate varies among site following a gamma distribution
#'       whose shape parameter is supplied via the gamma argument. By default
#'       (gamma = FALSE), no gamma correction is applied. Another assumption
#'       is that the base frequencies are balanced and thus equal to 0.25.
#'      \item
#'       K80: The distance derived by Kimura (1980), sometimes referred to as
#'       "Kimura's 2-parameters distance", has the same underlying assumptions
#'       than the Jukes-Cantor distance except that two kinds of substitutions
#'       are considered: transitions (A <-> G, C <-> T), and transversions
#'       (A <-> C, A <-> T, C <-> G, G <-> T). They are assumed to have
#'       different probabilities. A transition is the substitution of a purine
#'       (C, T) by another one, or the substitution of a pyrimidine (A, G) by
#'       another one. A transversion is the substitution of a purine by a
#'       pyrimidine, or vice-versa. Both transition and transversion rates are
#'       the same for all sites along the DNA sequence. Jin and Nei (1990)
#'       modified the Kimura model to allow for variation among sites following
#'       a gamma distribution. Like for the Jukes-Cantor model, the gamma
#'       shape parameter is supplied via the gamma argument. By default, no
#'       gamma correction is applied.
#'      \item
#'       F81: Felsenstein (1981) generalized the Jukes-Cantor model by relaxing
#'       the assumption of equal base frequencies. The formulae used in this
#'       function were taken from McGuire et al. (1999).
#'      \item
#'       K81: Kimura (1981) generalized his model (Kimura 1980) by assuming
#'       different rates for two kinds of transversions: A <-> C and G <-> T on
#'       one side, and A <-> T and C <-> G on the other. This is what Kimura
#'       called his "three substitution types model" (3ST), and is sometimes
#'       referred to as "Kimura's 3-parameters distance".
#'      \item
#'       F84: This model generalizes K80 by relaxing the assumption of equal
#'       base frequencies. It was first introduced by Felsenstein in 1984 in
#'       Phylip, and is fully described by Felsenstein and Churchill (1996).
#'       The formulae used in this function were taken from McGuire et al. (1999).
#'      \item
#'       BH87: Barry and Hartigan (1987) developed a distance based on the
#'       observed proportions of changes among the four bases. This distance
#'       is not symmetric. With by.pop = FALSE the result is returned as
#'       ape's full asymmetric matrix, not a dist object; with by.pop = TRUE
#'       the population averaging uses the block of distances from the
#'       individuals of the row population to those of the column population
#'       and assigns the mean to both triangles, silently symmetrising the
#'       result.
#'      \item
#'       T92: Tamura (1992) generalized the Kimura model by relaxing the
#'       assumption of equal base frequencies. This is done by taking into
#'       account the bias in G+C content in the sequences. The substitution
#'       rates are assumed to be the same for all sites along the DNA sequence.
#'      \item
#'       TN93: Tamura and Nei (1993) developed a model which assumes distinct
#'       rates for both kinds of transition (A <-> G versus C <-> T), and
#'       transversions. The base frequencies are not assumed to be equal and
#'       are estimated from the data. A gamma correction of the inter-site
#'       variation in substitution rates is possible, supplied via the gamma
#'       argument.
#'      \item
#'       GG95: Galtier and Gouy (1995) introduced a model where the G+C content
#'       may change through time. Different rates are assumed for transitons
#'       and transversions.
#'      \item
#'       logdet: The Log-Det distance, developed by Lockhart et al. (1994), is
#'       related to BH87. However, this distance is symmetric. Formulae from Gu
#'       and Li (1996) are used. dist.logdet in phangorn uses a different
#'       implementation that gives substantially different distances for
#'       low-diverging sequences.
#'      \item
#'       paralin: Lake (1994) developed the paralinear distance which can be
#'       viewed as another variant of the Barry-Hartigan distance.
#' }
#'
#' pairwise.missing : If TRUE, then missing values in the sequence (Ns and
#' ambiguity codes) will be accommodated in the calculations pair of taxa at
#' a time; otherwise, the deletion of data at positions in the sequence will
#' be global (deleted if any missing data at the position in any individual).
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param subst.model The evolutionary model of nucleotide substitutions to
#' employ in calculating genetic distance between individuals [default "F81"].
#' @param min.tag.len Minimum tag length of sequence tags to be used in the
#' analysis [default NULL].
#' @param pairwise.missing Whether to delete the sites with missing data in a
#' pairwise way [default TRUE].
#' @param by.pop If TRUE, the distance matrix is based on comparing
#' populations; if FALSE, on individuals [default TRUE].
#' @param gamma Gamma correction for inter-site variation in substitution
#' rates, passed to ape::dist.dna(): FALSE for no correction, or a positive
#' value giving the shape parameter (alpha) of the gamma distribution, for
#' the substitution models that support it [default FALSE].
#' @param variance If TRUE, ape::dist.dna() attaches the variances of the
#' individual-level distances as attribute 'variance'; the attribute is not
#' preserved by population averaging when by.pop = TRUE [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @return The distance matrix as an object of class dist. With
#' by.pop = TRUE, labels are popNames(x). With by.pop = FALSE, labels are
#' composites of the individual and population names in the form
#' indName_pop (e.g. 'T27_TENTERFIELD'), taken from the FASTA headers
#' written by gl2fasta(), not indNames(x). Exception: with
#' subst.model = "BH87" and by.pop = FALSE, the return is ape's full
#' asymmetric matrix (class matrix), not a dist object.
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#'
#' \donttest{
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' tmp <- gl.filter.monomorphs(platypus.gl, verbose = 0)
#' gl.dist.phylo(x=tmp,subst.model="F81")
#' }
#'
#' @import ape
#'
#' @export

gl.dist.phylo <- function(x,
                          subst.model = "F81",
                          min.tag.len = NULL,
                          pairwise.missing = TRUE,
                          by.pop = TRUE,
                          gamma = FALSE,
                          variance = FALSE,
                          verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v2025.1",
                   verbose = verbose)

  # STANDARD ERROR CHECKING
  # Check datatype
  datatype <- utils.check.datatype(x,
                                   accept = c("SNP", "SilicoDArT"),
                                   verbose = verbose)

  # Reject SilicoDArT before any work is done; sequence assembly requires
  # SNP data
  if (datatype == "SilicoDArT") {
    stop(
      error(
        "Fatal Error: Function gl.dist.phylo works only with SNP data. For SilicoDArT presence-absence data, use gl.dist.pop or gl.dist.ind\n"
      )
    )
  }

  if (!is(x, "dartR")) {
    class(x) <- "dartR"
    if (verbose > 2) {
      cat(
        warn(
          "Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"
        )
      )
      cat(
        warn(
          "                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR,
                 please use function dartR2gl\n"
        )
      )
    }
  }

  # FUNCTION SPECIFIC ERROR CHECKING

  # Population averaging needs at least two populations
  if (by.pop && nPop(x) < 2) {
    stop(
      error(
        "Fatal Error: by.pop = TRUE requires at least two populations to compare; this object has",
        nPop(x),
        "\n  Assign populations to the genlight object or set by.pop = FALSE for an individual-level matrix\n"
      )
    )
  }

  # Check for monomorphic loci and count heterozygous calls from the
  # genotype matrix directly -- gl.filter.monomorphs() errors on an
  # all-monomorphic object, and the sequence assembly densifies the
  # genotypes downstream in any case
  gm <- as.matrix(x)
  n.mono <- sum(apply(gm, 2, function(cc) {
    v <- cc[!is.na(cc)]
    length(unique(v)) <= 1 && !any(v == 1)
  }))
  if (n.mono == nLoc(x)) {
    stop(
      error(
        "Fatal Error: no polymorphic loci in the genlight object; distances cannot be computed. Provide data with polymorphic loci\n"
      )
    )
  }
  if (n.mono > 0) {
    if (verbose >= 2) {
      cat(warn("  Warning: genlight object contains", n.mono, "monomorphic loci\n"))
    }
  }

  # Heterozygous calls become IUPAC ambiguity codes in the assembled
  # sequences, which ape::dist.dna treats as missing data -- they carry no
  # distance signal. Warn at verbose >= 1 because the result is affected.
  n.called <- sum(!is.na(gm))
  n.het <- sum(gm == 1, na.rm = TRUE)
  if (n.het > 0 && verbose >= 1) {
    cat(
      warn(
        "  Warning:", n.het, "of", n.called, "genotype calls",
        paste0("(", round(100 * n.het / n.called, 1), "%)"),
        "are heterozygous and carry no distance signal\n"
      )
    )
    cat(
      warn(
        "    Heterozygotes are written as IUPAC ambiguity codes, which ape::dist.dna treats as missing data; distances rest on homozygous differences only\n"
      )
    )
    if (!pairwise.missing) {
      cat(
        warn(
          "    With pairwise.missing = FALSE, every site at which any individual is heterozygous (or missing) is deleted for ALL individuals, which can remove most or all variable sites\n"
        )
      )
    }
  }

  # DEFINE FUNCTIONS

  avg.dist <- function(gl, dist) {
    # A function to collapse a distance matrix to populations (by averaging)
    mat <- as.matrix(dist)
    # Number of populations or OTUs
    n_pop = nPop(gl)
    map <- pop(gl)
    names(map) <- indNames(gl)
    # Initialize an empty matrix to store aggregated distances
    aggregated_matrix = matrix(NA, nrow = n_pop, ncol = n_pop)
    rownames(aggregated_matrix) = colnames(aggregated_matrix) = popNames(gl)
    # Iterate over each pair of populations/OTUs
    for (i in 1:(n_pop - 1)) {
      for (j in (i + 1):n_pop) {
        # Indices of individuals in each population/OTU
        inds_i = which(map == rownames(aggregated_matrix)[i])
        inds_j = which(map == colnames(aggregated_matrix)[j])
        # Mean distance between populations/OTUs
        aggregated_matrix[i, j] = aggregated_matrix[j, i] = mean(mat[inds_i, inds_j], na.rm = TRUE)
      }
    }
    dist <- as.dist(aggregated_matrix)
    # dist <- as.matrix(dist)
    # rownames(dist) <- NULL
    # colnames(dist) <- NULL
    # names <- sprintf("%-10s", popNames(gl))
    # dist <- cbind(names,dist)
    return(dist)
  }

  # DO THE JOB

  if (!is.null(min.tag.len)) {
    if (verbose >= 2) {
      cat(report(
        "  Filtering sequence tags with length less than",
        min.tag.len,
        "\n"
      ))
    }
    x <- gl.filter.taglength(x, lower = min.tag.len, verbose = 0)
  }

  # Create the sequences in a form amenable to analysis by ape.
  # Full paths throughout -- no setwd, so a failure inside gl2fasta or
  # read.dna cannot strand the session working directory in tempdir()

  if (verbose >= 2) {
    cat(report("  Converting sequence tags to input format for package {ape}\n"))
  }

  gl2fasta(x,
           outfile = "tmp.fas",
           outpath = tempdir(),
           verbose = 0)

  sequences <- ape::read.dna(file.path(tempdir(), "tmp.fas"), format = "fasta")

  # Calculate distances between individuals
  if (verbose >= 2) {
    cat(report("  Calculating distances between individuals\n"))
    cat(report("    Substitution model:", subst.model, "\n"))
    cat(report("    Pairwise missing value deletion:", pairwise.missing, "\n"))
  }

  D <- ape::dist.dna(sequences,
                     model = subst.model,
                     gamma = gamma,
                     variance = variance,
                     pairwise.deletion = pairwise.missing)

  if (by.pop) {
    #Calculate average distances for pairwise populations
    if (verbose >= 2) {
      cat(report("  Calculating average distances between populations\n"))
      if (variance) {
        cat(warn("  Warning: the variance attribute applies to the individual-level distances and is not preserved by population averaging\n"))
      }
    }
    D <- avg.dist(gl = x, dist = D)

    #   hold <- getwd()
    #   setwd(tempdir())
    #     write(nPop(x),file="infile")
    #     write.table(D,file="infile",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
    #   setwd(hold)
    # } else {
    #   hold <- getwd()
    #   setwd(tempdir())
    #     write(nInd(x),file="infile")
    #     write.table(D,file="infile",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
    #   setwd(hold)
  }

  # Convert to a dist object
  # D <- D[,-1]
  # rownames(D) <- popNames(x)
  # colnames <- rownames
  # D <- as.dist(D)

  # FLAG SCRIPT END

  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }

  return(D)
}
