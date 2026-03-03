# # bootstrapping function
#' @name utils.allelic.richness
#' @title A utility script to calculate allelic richness using bootstraping
#' @family utilities
#' 
#' @param df Dataframe with SNP data [required].
#' @param indices indices [required].
#' @param boot_method Bootstrapping method
#' [default "loc"]
#'  
#' @author Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#'  
#' @export
#' @return calling function name

utils.allelic.richness <- function(df,
                                   indices,
                                   boot_method = "loc") {

  allelicRichness <- function(){}  #to hack package checking...
  
  Rcpp::cppFunction(
    '
double allelicRichness(NumericMatrix x) {
  int nrow = x.nrow();
  int ncol = x.ncol();

  // Vectors to store per-site values
  std::vector<int> rawCounts;      // total allele count per site
  std::vector<double> refTotals;   // sum of reference allele counts per site
  std::vector<double> altTotals;   // sum of alternate allele counts per site

  // Loop over sites (columns)
  for (int j = 0; j < ncol; j++) {
    int count0 = 0, count1 = 0, count2 = 0;
    // Loop over individuals (rows)
    for (int i = 0; i < nrow; i++) {
      double val = x(i, j);
      if (R_IsNA(val)) continue; // skip NA values
      if (val == 0) count0++;
      else if (val == 1) count1++;
      else if (val == 2) count2++;
    }
    int nonNA = count0 + count1 + count2;
    if(nonNA == 0) continue; // skip sites with no data
    int raw_count = 2 * nonNA; // each individual contributes 2 alleles
    double ref_total = 2.0 * count0 + 1.0 * count1;
    double alt_total = 2.0 * count2 + 1.0 * count1;

    rawCounts.push_back(raw_count);
    refTotals.push_back(ref_total);
    altTotals.push_back(alt_total);
  }

  int m = rawCounts.size();
  if (m == 0) return NA_REAL; // if no sites with data

  // Determine the minimum raw count (minimum sample size) across sites
  int min_sample_size = rawCounts[0];
  for (int i = 1; i < m; i++) {
    if (rawCounts[i] < min_sample_size)
      min_sample_size = rawCounts[i];
  }

  double sumRichness = 0.0;

  // For each site, compute r_ref and r_alt using the hypergeometric rarefaction formula:
  // r = 1 - choose(raw - allele_total, min_sample_size) / choose(raw, min_sample_size)
  // Instead of using choose() directly, we compute the ratio as a product.
  for (int i = 0; i < m; i++) {
    int r_count = rawCounts[i];
    double ref_total = refTotals[i];
    double alt_total = altTotals[i];

    double prod_ref = 1.0;
    // If (r_count - ref_total) is less than min_sample_size, then choose(...) is zero
    if ((r_count - ref_total) < min_sample_size) {
      prod_ref = 0.0;
    } else {
      for (int k = 0; k < min_sample_size; k++) {
        prod_ref *= ((r_count - ref_total - k) / static_cast<double>(r_count - k));
      }
    }
    double r_ref = 1.0 - prod_ref;

    double prod_alt = 1.0;
    if ((r_count - alt_total) < min_sample_size) {
      prod_alt = 0.0;
    } else {
      for (int k = 0; k < min_sample_size; k++) {
        prod_alt *= ((r_count - alt_total - k) / static_cast<double>(r_count - k));
      }
    }
    double r_alt = 1.0 - prod_alt;

    double site_richness = r_ref + r_alt;
    sumRichness += site_richness;
  }

  double meanRichness = sumRichness / m;
  // Round the result to 6 decimal places
  double factor = 1e6;
  double rounded = std::round(meanRichness * factor) / factor;

  return rounded;
}
    '
  )

  df <- df[indices,]

  if(boot_method == "loc"){
    df <- t(df)
  }

  res <- allelicRichness(df)

  return(res)

}
