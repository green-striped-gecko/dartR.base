% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.keep.pop.r
\name{gl.keep.pop}
\alias{gl.keep.pop}
\title{Removes all but the specified populations from a dartR genlight object}
\usage{
gl.keep.pop(
  x,
  pop.list,
  as.pop = NULL,
  recalc = FALSE,
  mono.rm = FALSE,
  verbose = NULL
)
}
\arguments{
\item{x}{Name of the genlight object [required].}

\item{pop.list}{List of populations to be retained [required].}

\item{as.pop}{Temporarily assign another locus metric as the population for
the purposes of deletions [default NULL].}

\item{recalc}{If TRUE, recalculate the locus metadata statistics [default FALSE].}

\item{mono.rm}{If TRUE, remove monomorphic and all NA loci [default FALSE].}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress but not results; 3, progress and results summary; 5, full report
 [default 2 or as specified using gl.set.verbosity].}
}
\value{
A reduced dartR genlight object
}
\description{
Individuals are assigned to populations based on associated specimen metadata
stored in the dartR genlight object. 
This script deletes all individuals apart from those in listed populations (pop.list).
Monomorphic loci and loci that are scored all NA are optionally deleted (mono.rm=TRUE). 
The script also optionally recalculates locus metatdata statistics to accommodate
the deletion of individuals from the dataset (recalc=TRUE).
The script returns a dartR genlight object with the retained populations 
and the recalculated locus metadata. The script works with both genlight objects
containing SNP genotypes and Tag P/A data (SilicoDArT).
}
\examples{
 # SNP data
   gl2 <- gl.keep.pop(testset.gl, pop.list=c('EmsubRopeMata', 'EmvicVictJasp'))
   gl2 <- gl.keep.pop(testset.gl, pop.list=c('EmsubRopeMata', 'EmvicVictJasp'),
   mono.rm=TRUE,recalc=TRUE)
   gl2 <- gl.keep.pop(testset.gl, pop.list=c('Female'),as.pop='sex')
 # Tag P/A data
   gs2 <- gl.keep.pop(testset.gs, pop.list=c('EmsubRopeMata','EmvicVictJasp'))
   
}
\seealso{
\code{\link{gl.drop.pop}} to drop rather than keep specified populations

Other data manipulation: 
\code{\link{gl.define.pop}()},
\code{\link{gl.drop.ind}()},
\code{\link{gl.drop.loc}()},
\code{\link{gl.drop.pop}()},
\code{\link{gl.edit.recode.pop}()},
\code{\link{gl.impute}()},
\code{\link{gl.join}()},
\code{\link{gl.keep.ind}()},
\code{\link{gl.keep.loc}()},
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
Custodian: Arthur Georges -- Post to
\url{https://groups.google.com/d/forum/dartr}
}
\concept{data manipulation}
