% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.make.recode.pop.r
\name{gl.make.recode.pop}
\alias{gl.make.recode.pop}
\title{Creates a proforma recode_pop_table file for reassigning population
 names
 @family data manipulation}
\usage{
gl.make.recode.pop(
  x,
  out.recode.file = "recode_pop_table.csv",
  outpath = NULL,
  verbose = NULL
)
}
\arguments{
\item{x}{Name of the genlight object [required].}

\item{out.recode.file}{File name of the output file (including extension)
[default recode_pop_table.csv].}

\item{outpath}{Directory to save the plot RDS files [default as specified 
by the global working directory or tempdir()]}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
[default 2 or as specified using gl.set.verbosity].}
}
\value{
A vector containing the new population names.
}
\description{
Renaming populations may be required when there have been errors in
assignment arising in the process from sample to sequence files or when one
wishes to amalgamate populations, or delete populations. Recoding populations
can also be done with a recode table (csv).
}
\details{
This function facilitates the construction of a recode table by producing a
proforma file with current population names in two identical columns. Edit
the second column to reassign populations. Use keyword 'Delete' to delete a
population.
When caution needs to be exercised because of the potential for breaking the
 'chain of evidence' associated with the samples, recoding individuals using
 a recode table (csv) can provide a clear record of the changes.
Use outpath=getwd() or when calling this function to direct output files 
to your working directory.
The function works with both genlight objects
containing SNP genotypes and Tag P/A data (SilicoDArT).
Apply the recoding using gl.recode.pop().
}
\examples{
result <- gl.make.recode.pop(testset.gl,out.recode.file='test.csv',outpath=tempdir(),verbose=2)
}
\author{
Custodian: Arthur Georges -- Post to
 \url{https://groups.google.com/d/forum/dartr}
}
