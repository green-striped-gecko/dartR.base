% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.report.diversity.r
\name{gl.report.diversity}
\alias{gl.report.diversity}
\title{Calculates diversity indexes for SNPs}
\usage{
gl.report.diversity(
  x,
  plot.display = TRUE,
  plot.theme = theme_dartR(),
  library = NULL,
  palette = NULL,
  plot.dir = NULL,
  plot.file = NULL,
  pbar = TRUE,
  table = "DH",
  verbose = NULL
)
}
\arguments{
\item{x}{Name of the genlight object containing the SNP or presence/absence
(SilicoDArT) data [required].}

\item{plot.display}{Specify if plot is to be displayed in the graphics window [default TRUE].}

\item{plot.theme}{User specified theme [default theme_dartR()].}

\item{library}{Name of the color library to be used [default scales::hue_pl].}

\item{palette}{Name of the color palette to be pulled from the specified library [default is library specific].}

\item{plot.dir}{Directory to save the plot RDS files [default as specified 
by the global working directory or tempdir()]}

\item{plot.file}{If TRUE, saves any ggplots and listings to the session
temporary directory (tempdir) [default FALSE].}

\item{pbar}{Report on progress. Silent if set to FALSE [default TRUE].}

\item{table}{Prints a tabular output to the console either 'D'=D values, or
'H'=H values or 'DH','HD'=both or 'N'=no table. [default 'DH'].}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
[default NULL, unless specified using gl.set.verbosity].}
}
\value{
A list of entropy indexes for each level of q and equivalent numbers
 for alpha and beta diversity.
}
\description{
This script takes a genlight object and calculates alpha and beta diversity
for q = 0:2. Formulas are taken from Sherwin et al. 2017. The paper describes
nicely the relationship between the different q levels and how they relate to
population genetic processes such as dispersal and selection. The citation 
below also includes a link to a 3-minute video that explains, q, D and H.
}
\details{
For all indexes, the entropies (H) and corresponding effective
 numbers, i.e. Hill numbers (D), which reflect the number of needed entities
 to get the observed values, are calculated. In a nutshell, the alpha indexes
 between the different q-values should be similar if there is no deviation
 from expected allele frequencies and occurrences (e.g. all loci in HWE &
 equilibrium). If there is a deviation of an index, this links to a process
 causing it, such as dispersal, selection or strong drift. For a detailed
 explanation of all the indexes, we recommend resorting to the literature
 provided below. Error bars are +/- 1 standard deviation.
 
\strong{ Function's output }

If the function's parameter "table" = "DH" (the default value) is used, the 
 output of the function is 20 tables.

The first two show the number of loci used. The name of each of the rest of 
the tables starts with three terms separated by underscores.

The first term refers to the q value (0 to 2).  The q values identify 
different ways of summarising diversity (H): q=0 is simply the number of 
alleles per locus, with no information about their relative proportions; q=2
 is the expected heterozygosity, ie the chance of drawing two different 
 alleles at random from the population; q=1 is the Shannon measure of 
 ‘surprise, relating to how likely it is that the next allele drawn will be
  one  that has not been seen before (Sherwin et al 2017, 2021, and 
  associated video). 

The second term refers to whether it is the diversity measure (H) or its 
transformation to Hill numbers (D)  The D value tells you how many 
equally-frequent alleles there would need to be to give the corresponding
 H-value (in the actual population)  The D-values are all in units of 
 numbers of alleles, so they can be plotted against the q-value to get a
  rich representation of the diversity  (Box 1, Fig II in Sherwin et al 
  2017, 2021, and associated video).

The third term refers to whether the diversity is calculated within 
populations (alpha) or between populations (beta). 

In the case of alpha diversity tables, standard deviations have their own 
table, which finishes with a fourth term: "sd".

In the case of beta diversity tables, standard deviations are in the upper 
triangle of the matrix and diversity values are in the lower triangle of the 
matrix.

\strong{ Plotting }

 Plot colours can be set with gl.select.colors().
 
 If plot.file is specified, plots are saved to the directory specified by the user, or the global
 default working directory set by gl.set.wd() or to the tempdir().

 Examples of other themes that can be used can be consulted in \itemize{
 \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
 \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
 }
}
\examples{
div <- gl.report.diversity(bandicoot.gl,library='brewer', table=FALSE,pbar=FALSE)
div$zero_H_alpha
div$two_H_beta
names(div)
}
\references{
Sherwin, W.B., Chao, A., Johst, L., Smouse, P.E. (2017, 2021). Information
 Theory Broadens the Spectrum of Molecular Ecology and Evolution. TREE
  32(12) 948-963. doi:10.1016/j.tree.2017.09.12 AND TREE 36:955-6 
  doi.org/10.1016/j.tree.2021.07.005 AND 3-Minute video:
  ars.els-cdn.com/content/image/1-s2.0-S0169534717302550-mmc2.mp4
}
\seealso{
Other unmatched report: 
\code{\link{gl.allele.freq}()},
\code{\link{gl.report.basics}()},
\code{\link{gl.report.heterozygosity}()}
}
\author{
Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr}),
 Contributors: William B. Sherwin, Alexander Sentinella
}
\concept{unmatched report}
