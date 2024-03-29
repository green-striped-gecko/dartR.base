% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl2bayescan.r
\name{gl2bayescan}
\alias{gl2bayescan}
\title{Converts a genlight object into a format suitable for input to Bayescan}
\usage{
gl2bayescan(x, outfile = "bayescan.txt", outpath = NULL, verbose = NULL)
}
\arguments{
\item{x}{Name of the genlight object containing the SNP data [required].}

\item{outfile}{File name of the output file (including extension)
[default bayescan.txt].}

\item{outpath}{Path where to save the output file [default global working 
directory or if not specified, tempdir()].}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
[default 2 or as specified using gl.set.verbosity].}
}
\value{
returns no value (i.e. NULL)
}
\description{
The output text file contains the SNP data and relevant BAyescan command
 lines to guide input.
}
\examples{
out <- gl2bayescan(testset.gl, outpath = tempdir())

}
\references{
Foll M and OE Gaggiotti (2008) A genome scan method to identify selected loci
 appropriate for both dominant and codominant markers: A Bayesian
  perspective. Genetics 180: 977-993.
}
\seealso{
Other linker: 
\code{\link{gl2bayesAss}()},
\code{\link{gl2bpp}()},
\code{\link{gl2demerelate}()},
\code{\link{gl2eigenstrat}()},
\code{\link{gl2faststructure}()},
\code{\link{gl2gds}()},
\code{\link{gl2genalex}()},
\code{\link{gl2genepop}()},
\code{\link{gl2geno}()},
\code{\link{gl2gi}()},
\code{\link{gl2hiphop}()},
\code{\link{gl2phylip}()},
\code{\link{gl2plink}()},
\code{\link{gl2related}()},
\code{\link{gl2sa}()},
\code{\link{gl2structure}()},
\code{\link{gl2treemix}()},
\code{\link{gl2vcf}()}
}
\author{
Custodian: Luis Mijangos (Post to
\url{https://groups.google.com/d/forum/dartr})
}
\concept{linker}
