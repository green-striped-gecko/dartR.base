% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.edit.recode.ind.r
\name{gl.edit.recode.ind}
\alias{gl.edit.recode.ind}
\title{Creates or edits individual (=specimen) names, creates a recode_ind
 file and applies the changes to a genlight object
 @family data manipulation}
\usage{
gl.edit.recode.ind(
  x,
  out.recode.file = NULL,
  outpath = NULL,
  recalc = FALSE,
  mono.rm = FALSE,
  verbose = NULL
)
}
\arguments{
\item{x}{Name of the genlight object [required].}

\item{out.recode.file}{Name of the file to output the new individual labels
[optional].}

\item{outpath}{Directory to save the plot RDS files [default as specified 
by the global working directory or tempdir()]}

\item{recalc}{If TRUE, recalculate the locus metadata statistics [default TRUE].}

\item{mono.rm}{If TRUE, remove monomorphic loci [default TRUE].}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress but not results; 3, progress and results summary; 5, full report
 [default 2 or as specified using gl.set.verbosity].}
}
\value{
An object of class ('genlight') with the revised individual labels.
}
\description{
A function to edit names of individual in a dartR genlight object, or to create a
reassignment table taking the individual labels from a genlight object, or to
edit existing individual labels in an existing recode_ind file. The amended 
recode table is then applied to the genlight object.
}
\details{
Renaming individuals may be required when there have been errors in labeling
 arising in the passage of samples to sequencing. There may be occasions
 where renaming individuals is required for preparation of figures. 
This function will input an existing recode table for editing and optionally
save it as a new table, or if the name of an input table is not supplied,
will generate a table using the individual labels in the parent genlight
object.
When caution needs to be exercised because of the potential for breaking the
'chain of evidence' associated with the samples, recoding individuals using
a recode table (csv) can provide a durable record of the changes.
For SNP genotype data, the function, having deleted individuals, optionally 
identifies resultant monomorphic loci or loci with all values missing 
and deletes them. The script also optionally recalculates the
locus metadata as appropriate. The optional deletion of monomorphic loci
and the optional recalculation of locus statistics is not available for
Tag P/A data (SilicoDArT).
Use outpath=getwd() when calling this function to direct
output files to your working directory.
The function returns a dartR genlight object with the new population assignments  
and the recalculated locus metadata.
}
\examples{
#this is an interactive example
if(interactive()){
gl <- gl.edit.recode.ind(testset.gl)
gl <- gl.edit.recode.ind(testset.gl, out.recode.file='ind.recode.table.csv')
}

}
\seealso{
\code{\link{gl.recode.ind}}, \code{\link{gl.drop.ind}},
\code{\link{gl.keep.ind}}
}
\author{
Custodian: Arthur Georges -- Post to
\url{https://groups.google.com/d/forum/dartr}
}
