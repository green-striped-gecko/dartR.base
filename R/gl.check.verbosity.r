#' @name gl.check.verbosity
#' @title Checks the current global verbosity
#' @family environment
#' 
#' @description
#' The verbosity can be set in one of two ways -- (a) explicitly by the user by
#' passing a value using the parameter verbose in a function, or (b) by setting
#' the verbosity globally as part of the r environment (gl.set.verbosity).

#' @param x User requested level of verbosity [default NULL].
#' 
#' @examples 
#' gl.check.verbosity()
#' 
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' 
#' @export
#' @return The verbosity, in variable verbose

# Version v.2023.3

# Function to check and set verbosity level
gl.check.verbosity <- function(x = NULL) {
  
  # SET VERBOSITY or GET it from global
  # If x is not provided, check if a global verbosity level has been set, otherwise use 2 as default.
  if (is.null(x)) {
    if (is.null(options()$dartR_verbose)) {
      verbose <- 2
    } else {
      verbose <- options()$dartR_verbose
    }
  } else {
    # If x is provided, check if it is a valid numeric value in the range of 0 to 5.
    if (is.numeric(x) & x >= 0 & x <= 5) {
      verbose <- x
    } else {
      # If x is not valid, issue a warning and use 2 as the default verbosity level.
      cat(warn("Warning: Parameter verbose must be an integer in the range 0 to 5, set to 2\n"))
      verbose <- 2
    }
  }
  
  # Return the determined verbosity level.
  return(verbose)
}

