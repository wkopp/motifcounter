\name{motif.length}
\alias{motif.length}
\title{Get the motif length}
\description{
This function returns length of the current motif. 
}
\usage{
motif.length()
}
\examples{

library(mdist)
pwmfile=system.file("extdata","x31.pwm", package="mdist")
read.motif(pwmfile,"tab")
motif.length()
}

