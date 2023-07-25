#' Setting up the package

#' Setting theme, colors and verbosity
#' @importFrom graphics axis barplot box image lines text
#' @importFrom grDevices hcl col2rgb rgb colorRampPalette rgb2hsv
#' @importFrom methods new is show getPackageName
#' @importFrom stats dist nobs optimize pchisq variable.names optim quantile rbinom cor var as.dist pnorm complete.cases
#' pgamma
#' @import ggplot2
#' @import dartR.data
#' @rawNamespace import(adegenet, except = plot)
#' @importFrom StAMPP stamppNeisD
#' @importFrom ape dist.gene
#' @import utils
#' @import patchwork
#' @import stringr
#' @importFrom crayon red yellow blue green cyan


zzz <- NULL  #to create a useful named help page

# defining function "dot" from data.table package to pass CRAN checks
`.` <- list

build = "v.2023.2"

# SET MESSAGES COLORS
# - For fatal errors use “error” which will print the message in red. Example
# usage: stop(error(“Fatal error”)) - For warning messages use “warn” which will
# print the message in yellow. Example usage: cat(warn(“message”)) - For
# reporting messages use “report” which will print the message in green. Example
# usage: cat(report(“message”)) - For important messages use “important” which
# will print the message in blue. Example usage: cat(important(“message”)) - For
# other messages as code use “code” which will print the message in cyan.
# Example usage: cat(code(“message”))
error <- crayon::red
warn <- crayon::yellow
report <- crayon::green
important <- crayon::blue
code <- crayon::cyan

# WELCOME MESSAGE 
.onAttach <- function(...) {
pn <- getPackageName()
packageStartupMessage(important(
  paste(
    "**** Welcome to",pn,"[Version",
    packageVersion(pn),
    "] ****\n"
  )
))
}
