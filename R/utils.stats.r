# Standard error
std.error <- function(x){
  sqrt(var(x, na.rm = TRUE) / sum(!is.na(x)))
  }
