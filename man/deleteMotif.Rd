\name{deleteMotif}
\alias{deleteMotif}
\title{Deletes the current motif}
\description{
This function unloads the current motif and frees the allocated memory 
}
\usage{
deleteMotif()
}
\examples{

library(mdist)
pwmfile=system.file("extdata","x31.pwm", package="mdist")
readMotif(pwmfile,"tab")
deleteMotif()
}

