% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.recalc.avgpic.r
\name{utils.recalc.avgpic}
\alias{utils.recalc.avgpic}
\title{A utility function to recalculate intermediate locus metrics}
\usage{
utils.recalc.avgpic(x, verbose = NULL)
}
\arguments{
\item{x}{Name of the genlight [required].}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
[default 2 or as specified using gl.set.verbosity]}
}
\value{
The modified genlight object.
}
\description{
WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.
}
\details{
Recalculates OneRatioRef, OneRatioSnp, PICRef, PICSnp,
 and AvgPIC by locus after some individuals or populations have been deleted.
 
The locus metadata supplied by DArT has OneRatioRef, OneRatioSnp, PICRef,
 PICSnp, and AvgPIC included, but the allelic composition will change when
 some individuals,or populations, are removed from the dataset and so the
 initial statistics will no longer apply. This script recalculates these
 statistics and places the recalculated values in the appropriate place in
 the genlight object.
If the locus metadata OneRatioRef|Snp, PICRef|Snp and/or AvgPIC do not exist,
 the script creates and populates them.
}
\seealso{
\code{utils.recalc.metrics} for recalculating all metrics,
\code{utils.recalc.callrate} for recalculating CallRate,
\code{utils.recalc.freqhomref} for recalculating frequency of homozygous
reference, \code{utils.recalc.freqhomsnp} for recalculating frequency of
homozygous alternate, \code{utils.recalc.freqhet} for recalculating frequency
 of heterozygotes, \code{gl.recalc.maf} for recalculating minor allele
 frequency, \code{gl.recalc.rdepth} for recalculating average read depth

Other utilities: 
\code{\link{gl.alf}()},
\code{\link{utils.check.datatype}()},
\code{\link{utils.dart2genlight}()},
\code{\link{utils.dist.binary}()},
\code{\link{utils.flag.start}()},
\code{\link{utils.hamming}()},
\code{\link{utils.het.pop}()},
\code{\link{utils.impute}},
\code{\link{utils.is.fixed}()},
\code{\link{utils.jackknife}()},
\code{\link{utils.n.var.invariant}()},
\code{\link{utils.plot.save}()},
\code{\link{utils.read.fasta}()},
\code{\link{utils.read.ped}()},
\code{\link{utils.recalc.callrate}()},
\code{\link{utils.recalc.freqhets}()},
\code{\link{utils.recalc.freqhomref}()},
\code{\link{utils.recalc.freqhomsnp}()},
\code{\link{utils.recalc.maf}()},
\code{\link{utils.reset.flags}()},
\code{\link{utils.transpose}()}
}
\author{
Custodian: Luis Mijangos (Post to
 \url{https://groups.google.com/d/forum/dartr})
}
\concept{utilities}
