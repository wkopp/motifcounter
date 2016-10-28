\name{motifLength}
\alias{motifLength}
\title{Get the motif length}
\description{
This function returns length of the current motif. 
}
\usage{
motifLength()
}
\examples{

library(mdist)
pwmfile=system.file("extdata", package="mdist")
readMotif(pwmfile)
motifLength()
}

