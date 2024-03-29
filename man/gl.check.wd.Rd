% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.check.wd.r
\name{gl.check.wd}
\alias{gl.check.wd}
\title{Checks the global working directory}
\usage{
gl.check.wd(wd = NULL, verbose = NULL)
}
\arguments{
\item{wd}{path to the working directory [default: tempdir()].}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
[default 2, unless specified using gl.set.verbosity].}
}
\value{
the working directory
}
\description{
The working directory can be set in one of two ways -- (a) explicitly by the user by
passing a value using the parameter plot.dir in a function, or (b) by setting
the working directory globally as part of the r environment (gl.setwd). The default is in acccordance to CRAN set to tempdir().
}
\examples{
gl.check.wd()

}
\seealso{
Other environment: 
\code{\link{gl.check.verbosity}()},
\code{\link{gl.print.history}()},
\code{\link{gl.set.wd}()},
\code{\link{theme_dartR}()}
}
\author{
Custodian: Bernd Gruber (Post to
\url{https://groups.google.com/d/forum/dartr})
}
\concept{environment}
