\name{delete.motif}
\alias{delete.motif}
\title{Deletes the current motif}
\description{
This function unloads the current motif. 
}
\usage{
delete.motif()
}
\examples{

library(mdist)
pwmfile=system.file("extdata","x31.pwm", package="mdist")
read.motif(pwmfile,"tab")
delete.motif()
}

