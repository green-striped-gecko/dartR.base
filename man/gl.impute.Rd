% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.impute.r
\name{gl.impute}
\alias{gl.impute}
\title{Imputes missing data}
\usage{
gl.impute(
  x,
  method = "neighbour",
  fill.residual = TRUE,
  parallel = FALSE,
  verbose = NULL
)
}
\arguments{
\item{x}{Name of the genlight object containing the SNP or presence-absence
data [required].}

\item{method}{Imputation method, either "frequency" or "HW" or "neighbour" 
or "random" [default "neighbour"].}

\item{fill.residual}{Should any residual missing values remaining after 
imputation be set to 0, 1, 2 at random, taking into account global allele 
frequencies at the particular locus [default TRUE].}

\item{parallel}{A logical indicating whether multiple cores -if available-
should be used for the computations (TRUE), or not (FALSE); requires the
package parallel to be installed [default FALSE].}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log ; 3, progress and results summary; 5, full report
[default 2 or as specified using gl.set.verbosity].}
}
\value{
A genlight object with the missing data imputed.
}
\description{
This function imputes genotypes on a population-by-population basis, where
populations can be considered panmictic, or imputes the state for
presence-absence data.
}
\details{
We recommend that imputation be performed on sampling locations, before
any aggregation. Imputation is achieved by replacing missing values using
either of four methods:
\itemize{
\item If "frequency", genotypes scored as missing at a locus in an individual
 are imputed using the average allele frequencies at that locus in the 
 population from which the individual was drawn.
\item If "HW", genotypes scored as missing at a locus in an individual are 
imputed by sampling at random assuming Hardy-Weinberg equilibrium. Applies 
only to genotype data.
\item If "neighbour", substitute the missing values for the focal individual
 with the values taken from the nearest neighbour. Repeat with next nearest
 and so on until all missing values are replaced.
\item if "random", missing data are substituted by random values (0, 1 or 2). 
}
  The nearest neighbour is the one at the smallest Euclidean distancefrom 
  the focal individual
  The advantage of this approach is that it works regardless of how many
  individuals are in the population to which the focal individual belongs,
  and the displacement of the individual is haphazard as opposed to:
  (a) Drawing the individual toward the population centroid (HW and Frequency).
  (b) Drawing the individual toward the global centroid (glPCA).
Note that loci that are missing for all individuals in a population are not 
imputed with method 'frequency' or 'HW' and can give unpredictable results
for particular individuals using 'neighbour'. Consider using the function 
\code{\link{gl.filter.allna}} with by.pop=TRUE to remove them first.
}
\examples{
 \donttest{
require("dartR.data")
# SNP genotype data
gl <- gl.filter.callrate(platypus.gl,threshold=0.95)
gl <- gl.filter.allna(gl)
gl <- gl.impute(gl,method="neighbour")
# Sequence Tag presence-absence data
gs <- gl.filter.callrate(testset.gs,threshold=0.95)
gl <- gl.filter.allna(gl)
gs <- gl.impute(gs, method="neighbour")
}
gs <- gl.impute(platypus.gl,method ="random")

}
\seealso{
Other data manipulation: 
\code{\link{gl.define.pop}()},
\code{\link{gl.drop.ind}()},
\code{\link{gl.drop.loc}()},
\code{\link{gl.drop.pop}()},
\code{\link{gl.edit.recode.pop}()},
\code{\link{gl.join}()},
\code{\link{gl.keep.ind}()},
\code{\link{gl.keep.loc}()},
\code{\link{gl.keep.pop}()},
\code{\link{gl.make.recode.ind}()},
\code{\link{gl.merge.pop}()},
\code{\link{gl.reassign.pop}()},
\code{\link{gl.recode.ind}()},
\code{\link{gl.recode.pop}()},
\code{\link{gl.rename.pop}()},
\code{\link{gl.sample}()},
\code{\link{gl.sim.genotypes}()},
\code{\link{gl.sort}()},
\code{\link{gl.subsample.ind}()},
\code{\link{gl.subsample.loc}()}
}
\author{
Custodian: Luis Mijangos 
(Post to \url{https://groups.google.com/d/forum/dartr})
}
\concept{data manipulation}
