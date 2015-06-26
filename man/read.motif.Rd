\name{read.motif}
\alias{read.motif}
\title{Read and load the PWM}
\description{
The function reads the PWM from a file. The PWM can be encoded
in various formats.
}
\usage{
read.motif(pwmfile,format, pseudocount=0.01)
}
\arguments{
\item{pwmfile}{Character string of the name of the file which contains the PWM. }

\item{format}{Character string which determines the format of the file. 
Currently "tab" and "transfac" are supported.}
\item{pseudocount}{Numeric value which is added to each element of the loaded 
PWM in order to avoid zero probability elements. The columns of the matrix are
renormalized automatically.}
}
\examples{

library(mdist)
pwmfile=system.file("extdata","x31.pwm", package="mdist")
read.motif(pwmfile,"tab", 0.01)
}

