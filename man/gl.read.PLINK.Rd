% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.read.PLINK.r
\name{gl.read.PLINK}
\alias{gl.read.PLINK}
\title{Reads PLINK data file into a genlight object}
\usage{
gl.read.PLINK(
  filename,
  ind.metafile = NULL,
  loc.metafile = NULL,
  plink.cmd = "plink",
  plink.path = "path",
  plink.flags = NULL,
  verbose = NULL
)
}
\arguments{
\item{filename}{Fully qualified path to PLINK input file (without including the extension)}

\item{ind.metafile}{Name of the csv file containing the metrics for
individuals [optional].}

\item{loc.metafile}{Name of the csv file containing the metrics for
loci [optional].}

\item{plink.cmd}{The 'name' to call plink. This will depend on the file name 
(without the extension '.exe' if on windows) or the name of the PATH variable}

\item{plink.path}{The path where the executable is. If plink is listed in
the PATH then there is no need for this. This is what the option "path"
means}

\item{plink.flags}{additional possible parameters passed on to plink.}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
[default 2 or as specified using gl.set.verbosity].}
}
\value{
A genlight object with the SNP data and associated metadata included.
}
\description{
This function imports PLINK data into a genlight object and append 
available metadata.
}
\details{
This function handles .ped or .bed file (with the associate files - 
e.g. .fam, .bim). However, if a .ped file is provided, PLINK needs to be 
installed and it is used to convert the .ped into a .bed, which is then
 converted into a genlight.

Additional metadata can be included passing .csv files. These will be appended 
to the existing metadata present in the PLINK files. 

The locus metadata needs to be in a csv file with headings, with a mandatory
column headed AlleleID corresponding exactly to the locus identity labels
provided with the SNP data
}
\author{
Custodian: Carlo Pacioni -- Post to
\url{https://groups.google.com/d/forum/dartr}
}
