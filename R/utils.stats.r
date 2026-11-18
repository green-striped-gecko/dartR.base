#' Standard error of the mean: sqrt(var/n) over the non-missing values
#' @param x Numeric vector (NAs excluded)
#' @return The standard error of the mean
#' @noRd
std.error <- function(x){
  sqrt(var(x, na.rm = TRUE) / sum(!is.na(x)))
  }
