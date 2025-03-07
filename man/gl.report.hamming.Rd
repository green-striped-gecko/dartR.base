% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.report.hamming.r
\name{gl.report.hamming}
\alias{gl.report.hamming}
\title{Calculates the pairwise Hamming distance between DArT trimmed DNA
sequences}
\usage{
gl.report.hamming(
  x,
  rs = 5,
  threshold = 3,
  tag.length = 69,
  plot.display = TRUE,
  plot.theme = theme_dartR(),
  plot.colors = NULL,
  plot.dir = NULL,
  plot.file = NULL,
  probar = FALSE,
  verbose = NULL
)
}
\arguments{
\item{x}{Name of the genlight object containing the SNP data [required].}

\item{rs}{Number of bases in the restriction enzyme recognition sequence
[default 5].}

\item{threshold}{Minimum acceptable base pair difference for display on the
boxplot and histogram [default 3].}

\item{tag.length}{Typical length of the sequence tags [default 69].}

\item{plot.display}{Specify if plot is to be produced [default TRUE].}

\item{plot.theme}{User specified theme [default theme_dartR()].}

\item{plot.colors}{Vector with two color names for the borders and fill
[default c("#2171B5", "#6BAED6")].}

\item{plot.dir}{Directory to save the plot RDS files [default as specified 
by the global working directory or tempdir()]}

\item{plot.file}{Filename (minus extension) for the RDS plot file [Required for plot save]}

\item{probar}{If TRUE, a progress bar is displayed during run [defalut FALSE]}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
[default 2, unless specified using gl.set.verbosity].}
}
\value{
Returns unaltered genlight object
}
\description{
Hamming distance is calculated as the number of base differences
between two sequences which can be expressed as a count or a proportion.
Typically, it is calculated between two sequences of equal length. In the
context of DArT trimmed sequences, which differ in length but which are
anchored to the left by the restriction enzyme recognition sequence, it is
sensible to compare the two trimmed sequences starting from immediately after
the common recognition sequence and terminating at the last base of the
shorter sequence.
}
\details{
The function \code{\link{gl.filter.hamming}} will filter out one of
two loci if their Hamming distance is less than a specified percentage
Hamming distance can be computed by exploiting the fact that the dot product
of two binary vectors x and (1-y) counts the corresponding elements that are
different between x and y. This approach can also be used for vectors that
contain more than two possible values at each position (e.g. A, C, T or G).
If a pair of DNA sequences are of differing length, the longer is truncated.
The algorithm is that of Johann de Jong
\url{https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/}
as implemented in \code{\link{utils.hamming}}
  If plot.file is specified, plots are saved to the directory specified by the user, or the global
  default working directory set by gl.set.wd() or to the tempdir().
 Examples of other themes that can be used can be consulted in \itemize{
 \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
 \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
 }
}
\examples{
 \donttest{
gl.report.hamming(testset.gl[,1:100])
gl.report.hamming(testset.gs[,1:100])
}
#' # SNP data
test <- platypus.gl
test <- gl.subsample.loc(platypus.gl,n=50)
result <- gl.report.hamming(test, verbose=3)
result <- gl.report.hamming(test, plot.file="ttest", verbose=3)
}
\seealso{
\code{\link{gl.filter.hamming}}

Other matched report: 
\code{\link{gl.filter.excess.het}()},
\code{\link{gl.report.allna}()},
\code{\link{gl.report.callrate}()},
\code{\link{gl.report.locmetric}()},
\code{\link{gl.report.maf}()},
\code{\link{gl.report.overshoot}()},
\code{\link{gl.report.pa}()},
\code{\link{gl.report.rdepth}()},
\code{\link{gl.report.reproducibility}()},
\code{\link{gl.report.secondaries}()},
\code{\link{gl.report.taglength}()}
}
\author{
Custodian: Arthur Georges -- Post to 
\url{https://groups.google.com/d/forum/dartr}
}
\concept{matched report}
