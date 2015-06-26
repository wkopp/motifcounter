\name{matrix2motif}
\alias{matrix2motif}
\title{Loads the PWM model from an R matrix}
\description{
This function loads the PWM model of interest into the memory.
}
\usage{
matrix2motif(data)
}
\arguments{
\item{data}{Must be an 4xM PWM matrix where M is the length of the motif.
The columns of the matrix must be normalized so that they correspond to
probabilities.}

}
\examples{

library(mdist)
pwmfile=system.file("extdata","x31.pwm", package="mdist")
m=t(as.matrix(read.table(pwmfile)))
#read.motif(pwmfile,"tab")
matrix2motif(m)
}


