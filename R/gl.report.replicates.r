#' @name gl.report.replicates
#' @title Identify replicated individuals 
#' @description
#' This function scans a genlight object for pairs of individuals that
#' share a high proportion of identical genotype calls across a minimum number of
#' jointly observed (non-missing) loci. It is intended to help detect technical
#' replicates, sample duplicates, and other near-identical samples, and to suggest
#' which individual to remove based on missing-data rate.
#' @param x Name of the genlight object containing the SNP data [required].
#' @param loc_threshold Minimum number of loci required to asses that two 
#' individuals are replicates [default 100].
#' @param perc_geno Minimum percentage of genotypes in which two individuals 
#' should be the same [default 0.95]. 
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot_theme User specified theme [default theme_dartR()].
#' @param plot_colors Vector with two color names for the borders and fill
#' [default c("#2171B5", "#6BAED6")].
#' @param bins Number of bins to display in histograms [default 100].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' ## What the function computes
#' 
#' The input is coerced to a numeric matrix. For every pair of
#' individuals \eqn{i,j}, the function counts:
#' \itemize{
#'   \item `nloc` = the number of loci where both individuals have non-missing
#'     genotype calls (pairwise complete observations).
#'   \item `nsame` = the number of those `nloc` loci where the two genotype calls
#'     are identical.
#'   \item `perc` = `nsame / nloc`, the proportion of identical genotype calls for
#'     that pair.
#' }
#' Pairs are reported as putative replicates/duplicates when both conditions hold:
#' `nloc > loc_threshold` and `perc > perc_geno`.
#'
#' ## How the computation is implemented 
#' 
#' The pairwise counts (`nsame` and `nloc`) are computed with a C++ routine 
#' compiled at run time. This requires the Rcpp/RcppParallel packages. The 
#' first call in a session will be slower because the C++ code must be compiled;
#' subsequent calls are fast.
#'
#' ## Suggested individual to drop
#' For each flagged pair, the function calculates per-individual missingness as the
#' proportion of `NA` genotypes across all loci. The replicate suggested for removal
#' (`ind_to_drop`) is the member of the pair with the higher missingness, so that
#' the retained replicate maximises usable genotype information.
#'
#' ## How to interpret the histograms
#' When `plot.out = TRUE`, the function draws:
#' \enumerate{
#'   \item A histogram of `perc` across *all* pairwise comparisons.
#'   \item A zoomed histogram restricted to `perc > 0.8` to resolve the upper tail.
#' }
#' In an ideal large dataset containing a mixture of unrelated and related
#' individuals plus several replicated individuals (e.g., capture/mark/recapture),
#' the first histogram often shows four approximate modes (“peaks”):
#' \enumerate{
#'   \item Unrelated pairs (lowest `perc`).
#'   \item Second-degree relatives (e.g., cousins).
#'   \item First-degree relatives (e.g., parent–offspring, full siblings).
#'   \item Replicated/duplicate samples (highest `perc`, near 1).
#' }
#' To confidently separate true replicates from close relatives, the zoomed
#' histogram should show a clear gap between the third and fourth peaks: one or
#' more bins with zero counts between the close-relative peak and the replicate
#' peak. If the close-relative and replicate distributions overlap (no zero-count
#' gap), you should treat replicate calls as uncertain and consider increasing
#' marker quality/quantity, tightening filters, or raising `loc_threshold` and/or
#' `perc_geno`.
#' @return A list with three elements:
#'\itemize{
#'\item table.rep: A dataframe with pairwise results of percentage of same 
#'genotypes between two individuals, the number of loci used in the comparison 
#'and the missing data for each individual.
#'\item ind.list.drop: A vector of replicated individuals to be dropped.
#' Replicated individual with the least missing data is reported.
#'\item ind.list.rep: A list of of each individual that has replicates in the 
#'dataset, the name of the replicates and the percentage of the same genotype.
#'  }
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \donttest{
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' res_rep <- gl.report.replicates(platypus.gl, loc_threshold = 500, 
#' perc_geno = 0.85)
#' }
#' @family report functions
#' @export

gl.report.replicates <- function(x,
                                 loc_threshold = 100,    
                                 perc_geno = 0.95,       
                                 plot.out = TRUE,        
                                 plot_theme = theme_dartR(), 
                                 plot_colors = c("#2171B5", "#6BAED6"),
                                 bins = 100,             
                                 verbose = NULL){
  
  # Pre-allocate variables to avoid "no visible binding" notes
  ind1 <- ind1_miss <- ind2 <- ind2_miss <- ind_to_drop <- nloc <- perc <- NULL
  
  # Determine verbosity level (internal helper)
  verbose  <- gl.check.verbosity(verbose)
  
  # Record function name for logging
  funname  <- match.call()[[1]]
  
  # Flag the start of this function call (internal helper)
  utils.flag.start(func = funname, build = "Jody", verbose = verbose)
  
  # Convert genlight or other object to plain numeric matrix
  xx <- as.matrix(x)
  
  # Stub to satisfy package checks when defining the real function below
  pairwiseMatchParallel <- function() {}
  
  # Define and compile a C++ function that computes pairwise shared/non-missing 
  # counts in parallel
  suppressWarnings(
  Rcpp::cppFunction(
    code = '
    #include <Rcpp.h>
    #include <RcppParallel.h>
    using namespace Rcpp;
    using namespace RcppParallel;

    // Worker struct for computing pairwise matches in parallel
    struct PairwiseWorker : public Worker {
      const RMatrix<double> x;
      RMatrix<double>       same;     // counts of identical alleles
      RMatrix<double>       nonmiss;  // counts of non-missing comparisons
      const int             m;        // number of loci (columns)

      // Constructor: bind input matrix and result matrices
      PairwiseWorker(const NumericMatrix& x_,
                     NumericMatrix&       same_,
                     NumericMatrix&       nonmiss_)
        : x(x_), same(same_), nonmiss(nonmiss_), m(x_.ncol()) {}

      // The work to do on each chunk of rows
      void operator()(std::size_t begin, std::size_t end) override {
        std::size_t n = x.nrow();
        for (std::size_t i = begin; i < end; ++i) {
          for (std::size_t j = i + 1; j < n; ++j) {
            int nsame = 0, nobs = 0;
            // Loop over loci for this pair
            for (int k = 0; k < m; ++k) {
              double xi = x(i,k), xj = x(j,k);
              // Count only non-missing comparisons
              if (!traits::is_na<REALSXP>(xi) &&
                  !traits::is_na<REALSXP>(xj)) {
                ++nobs;
                if (xi == xj) ++nsame;
              }
            }
            // Fill both [i,j] and [j,i] entries
            same(i,j)    = same(j,i)    = nsame;
            nonmiss(i,j) = nonmiss(j,i) = nobs;
          }
        }
      }
    };

    // [[Rcpp::export]]
    List pairwiseMatchParallel(const NumericMatrix& x) {
      int n = x.nrow();
      NumericMatrix same(n,n), nonmiss(n,n);
      PairwiseWorker worker(x, same, nonmiss);
      // Parallel loop over rows 0 to n-1
      parallelFor(0, n-1, worker);
      return List::create(_["same"]=same, _["nonmiss"]=nonmiss);
    }
  ',
    depends = "RcppParallel",
    plugin  = "cpp11"
  )
  )
  
  # Call the compiled C++ function
  pm <- pairwiseMatchParallel(xx)
  mat_same     <- pm[["same"]]      # matrix of identical counts
  mat_nonmiss  <- pm[["nonmiss"]]   # matrix of non-missing counts
  
  # Assign individual names to rows/columns
  dimnames(mat_same)    <- dimnames(mat_nonmiss) <-
    list(indNames(x), indNames(x))
  
  # Compute proportion identical at each pair
  mat_prop <- mat_same / mat_nonmiss
  
  # Melt into a long table: one row per pair
  tab <- data.table(
    ind1 = rep(indNames(x), each = nrow(mat_prop)),
    ind2 = rep(indNames(x),  times = nrow(mat_prop)),
    perc = as.vector(mat_prop),     # proportion identical
    nloc = as.vector(mat_nonmiss)   # number of loci compared
  )[!is.na(perc)]                   # drop any NA proportions
  
  # Find pairs exceeding both thresholds
  col_same <- tab[nloc > loc_threshold & perc > perc_geno][order(-perc)]
  
  # If none, return a message suggesting to lower thresholds
  if (!nrow(col_same)) {
    msg <- sprintf(
      "No pair of individuals share > %.1f%% identical genotypes across > %d loci.  Lower the thresholds?",
      perc_geno * 100, loc_threshold
    )
    return(msg)
  }
  
  # Calculate per-individual missing-data proportion
  miss_prop <- 1 - rowSums(!is.na(xx)) / ncol(xx)
  miss_dt   <- data.table(id = indNames(x), miss = miss_prop)
  
  # Merge in missing-data rates for ind1 and ind2
  col_same <- merge(col_same, miss_dt, by.x = "ind1", by.y = "id")
  setnames(col_same, "miss", "ind1_miss")
  col_same <- merge(col_same, miss_dt, by.x = "ind2", by.y = "id")
  setnames(col_same, "miss", "ind2_miss")
  
  # Decide which replicate to drop: the one with higher missing proportion
  col_same[, ind_to_drop := ifelse(ind1_miss > ind2_miss, ind1, ind2)]
  ind_list <- unique(col_same$ind_to_drop)
  
  # Also prepare a list of all replicates per focal sample
  keep_pairs   <- mat_prop >= perc_geno & mat_nonmiss > loc_threshold
  ind_list_rep <- apply(keep_pairs, 2, function(v) indNames(x)[v])
  ind_list_rep <- ind_list_rep[lengths(ind_list_rep) > 0]
  
  # If plotting is requested, show overall and zoomed-in histograms
  if (plot.out) {
    p_all <- ggplot(tab, aes(perc)) +
      geom_histogram(bins = bins,
                     colour = plot_colors[1],
                     fill   = plot_colors[2]) +
      ylab("Count") + xlab("") + plot_theme
    
    p_zoom <- ggplot(tab[perc > 0.8], aes(perc)) +
      geom_histogram(bins = bins,
                     colour = plot_colors[1],
                     fill   = plot_colors[2]) +
      ylab("Count") + xlab("Proportion identical") + plot_theme
    
    # Print the two plots stacked vertically
    print(p_all / p_zoom)
  }
  
  # Final verbose message if requested
  if (verbose >= 1) cat(report("Completed:", funname, "\n"))
  
  list(
    table.rep     = col_same[],
    ind.list.drop = ind_list,
    ind.list.rep  = ind_list_rep
  )
  
}
