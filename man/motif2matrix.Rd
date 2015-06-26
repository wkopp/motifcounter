\name{motif2matrix}
\alias{motif2matrix}
\title{Get the motif in an R matrix}
\description{
This function returns the current PFM as an R matrix.
}
\usage{
motif2matrix()
}
\examples{

library(mdist)
pwmfile=system.file("extdata","x31.pwm", package="mdist")
read.motif(pwmfile,"tab")
x=motif2matrix()
}

