\name{mdist.option}
\alias{mdist.option}
\title{Set parameters for the enrichment analysis}
\description{
This function sets parameters for the enrichment analysis.
}
\usage{
mdist.option(alpha, gran, ncores=1)
}
\arguments{
\item{alpha}{Numerical value that defines the significance level for
calling motif hits.}

\item{gran}{Numerical value that defines the granularity which shall be
used for discretizing the score  range.}
\item{ncores}{Number of cores used for parallel processing. Default: ncores=1.}
}
\examples{

library(mdist)
mdist.option(0.01,0.01)
}

