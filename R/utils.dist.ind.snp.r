#' @name utils.dist.ind.snp
#' @title Calculates a distance matrix for individuals defined in a
#' genlight object using SNP data (DArTseq)
#' @family distance
#'
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' @param x Name of the genlight containing the genotypes [required].
#' @param method Specify distance measure [default Euclidean].
#' @param scale If TRUE and method='Euclidean', the summed squared
#' differences are divided by the number of loci compared, so the distance
#' falls in the range [0,2] for SNP genotypes scored 0, 1 or 2
#' [default FALSE].
#' @param type Specify the format and class of the object to be returned,
#' dist for a object of class dist, matrix for an object of class matrix
#'  [default "dist"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log; 3, progress and results summary; 5, full report [default 2].
#'
#' @details
#' This script calculates various distances between individuals based on
#' SNP genotypes.
#'
#' The distance measure can be one of:
#'  \itemize{
#'   \item Euclidean -- Euclidean Distance applied to Cartesian coordinates
#'    defined by the loci, scored as 0, 1 or 2.
#'  \item Simple -- based on shared alleles per locus: 0 where no alleles
#'  are shared, 1 where one allele is shared, 2 where both alleles are
#'  shared; distance = 1 - mean(shared)/2. With biallelic SNPs this is
#'  arithmetically identical to the Czekanowski distance.
#'  \item Absolute -- 0 where no alleles are shared (opposite homozygotes),
#'  1 where one or both alleles are shared; distance is the proportion of
#'  loci with no shared alleles.
#'  \item Czekanowski (or Manhattan) calculates the city block metric distance
#'  by summing the scores on each axis (locus).
#'  }
#'
#' All measures are invariant to which allele is scored as the reference
#' allele.
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @keywords internal
# @export
#' @return An object of class 'dist' or 'matrix' giving distances between
#'  individuals

# Examples for testing
# D <- utils.dist.ind.snp(testset.gl, method='Manhattan')
# D <- utils.dist.ind.snp(testset.gl, method='Euclidean', scale = TRUE)
# D <- utils.dist.ind.snp(testset.gl, method='Simple')

utils.dist.ind.snp <- function(x,
                               method  = "Euclidean",
                               scale   = FALSE,
                               type    = "dist",
                               verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  funname <- utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION-SPECIFIC SETTINGS
  method <- tolower(method)

  if (!(method %in% c("euclidean", "simple", "absolute", "czekanowski", "manhattan"))) {
    if (verbose >= 2) {
      cat(warn("  Warning: Method not in the list of options, set to Euclidean\n"))
    }
    method <- "euclidean"
  }

  if (scale && method != "euclidean") {
    if (verbose >= 1) {
      cat(warn("  Warning: parameter scale only applies to Euclidean Distance, ignored\n"))
    }
  }

  # DO THE JOB
  if (verbose >= 2) {
    if (method == "euclidean") {
      if (scale) {
        cat(report("  Calculating the scaled distance matrix --", method, "\n"))
      } else {
        cat(report("  Calculating the unscaled distance matrix --", method, "\n"))
      }
    } else {
      cat(report("  Calculating the distance matrix --", method, "\n"))
    }
  }

  if (.has_fbm(x)) {
    # File-backed: the per-pair counts are added over blocks of loci read
    # from the backing file (.gl_pair_counts), so the genotype matrix is
    # never decoded as a whole. The counts are integers and the distances
    # are the same expressions as in dist_mod(), so the result is
    # identical to that of the gen-backed object.
    dd <- .dist_ind_snp_fbm(x, method = method, scale = scale)
  } else {
  mat <- as.matrix(x)

  # dummy so R CMD check does not flag the Rcpp-compiled symbol as
  # undefined; the real function is created by cppFunction below
  dist_mod <- function() {}
  Rcpp::cppFunction(
'
    // Compute a symmetric distance matrix for rows of x
    // method: euclidean, simple, absolute, manhattan or czekanowski
    // scale: if true, normalize Euclidean distances by number of loci
    NumericMatrix dist_mod(const NumericMatrix& x,
                           std::string method = "euclidean",
                           bool scale = true) {
      // All measures are symmetric in the allele labelling: they depend
      // on the genotypes only through the number of unshared alleles, so
      // recoding the reference allele leaves them unchanged.
      int nI = x.nrow();
      int nL = x.ncol();

      NumericMatrix dd(nI, nI);

      for (int i = 0; i < nI; ++i) {
        dd(i, i) = 0.0;
        for (int j = i + 1; j < nI; ++j) {
          double val = 0.0;

          if (method == "euclidean") {
            double sumsq = 0.0;
            int L = 0;
            for (int k = 0; k < nL; ++k) {
              double v1 = x(i, k), v2 = x(j, k);
              if (NumericMatrix::is_na(v1) || NumericMatrix::is_na(v2)) continue;
              double d = v1 - v2;
              sumsq += d * d;
              ++L;
            }
            if (L > 0) {
              val = scale ? std::sqrt(sumsq / L) : std::sqrt(sumsq);
            } else {
              val = NA_REAL;
            }

          } else if (method == "simple") {
            // shared alleles per locus: S = 2 - |g1 - g2| (0, 1 or 2);
            // distance = 1 - mean(S)/2
            double sumshared = 0.0;
            int L = 0;
            for (int k = 0; k < nL; ++k) {
              double v1 = x(i, k), v2 = x(j, k);
              if (NumericMatrix::is_na(v1) || NumericMatrix::is_na(v2)) continue;
              sumshared += 2.0 - std::abs(v1 - v2);
              ++L;
            }
            val = (L > 0 ? 1.0 - sumshared / (2.0 * L) : NA_REAL);

          } else if (method == "absolute") {
            // shared-any per locus: 1 unless opposite homozygotes;
            // distance = proportion of loci sharing no alleles
            double sumany = 0.0;
            int L = 0;
            for (int k = 0; k < nL; ++k) {
              double v1 = x(i, k), v2 = x(j, k);
              if (NumericMatrix::is_na(v1) || NumericMatrix::is_na(v2)) continue;
              sumany += (std::abs(v1 - v2) < 2.0 ? 1.0 : 0.0);
              ++L;
            }
            val = (L > 0 ? 1.0 - sumany / L : NA_REAL);

          } else if (method == "manhattan" || method == "czekanowski") {
            double sumabs = 0.0;
            int L = 0;
            for (int k = 0; k < nL; ++k) {
              double v1 = x(i, k), v2 = x(j, k);
              if (NumericMatrix::is_na(v1) || NumericMatrix::is_na(v2)) continue;
              sumabs += std::abs(v1 - v2);
              ++L;
            }
            val = (L > 0 ? sumabs / (2.0 * L) : NA_REAL);

          } else {
            val = NA_REAL;
          }

          dd(i, j) = val;
          dd(j, i) = val;
        }
      }

      return dd;
    }'
  )

  dd <- dist_mod(mat, method = method, scale = scale)
  }

  if (type == "dist") {
    dd <- as.dist(dd)
    if (verbose >= 2) {
      cat(report("  Returning a stats::dist object\n"))
    }
  } else {
    if (verbose >= 2) {
      cat(report("  Returning a square matrix object\n"))
    }
  }

  # FLAG SCRIPT END
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }

  return(dd)
}

## Distances of utils.dist.ind.snp() for a file-backed object, from the
## pairwise counts of .gl_pair_counts(): L loci scored in both individuals,
## n1 differing by one and n2 by two allele copies. With |g1 - g2| in
## {0, 1, 2}: sum of squared differences n1 + 4 n2, sum of absolute
## differences n1 + 2 n2, shared alleles 2 L - (n1 + 2 n2), loci sharing an
## allele L - n2. All sums are whole numbers, as in dist_mod(); a pair
## without common loci is NA, the diagonal 0.
.dist_ind_snp_fbm <- function(x, method, scale) {
  pc <- .gl_pair_counts(x)
  L <- pc$L
  sumabs <- pc$n1 + 2 * pc$n2
  v <- switch(method,
    euclidean = {
      sumsq <- pc$n1 + 4 * pc$n2
      if (scale) sqrt(sumsq / L) else sqrt(sumsq)
    },
    simple = 1 - (2 * L - sumabs) / (2 * L),
    absolute = 1 - (L - pc$n2) / L,
    manhattan = ,
    czekanowski = sumabs / (2 * L)
  )
  v[L == 0] <- NA_real_
  .gl_dist_to_matrix(v, nInd(x))
}
