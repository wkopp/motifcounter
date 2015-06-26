\name{store.background}
\alias{store.background}
\title{Store background model for later use in a binary file}
\description{
This function stores the background model, in order for it to be reused
later.
}
\usage{
store.background(file)
}
\arguments{
\item{file}{Name of the binary file, where the background model should
be stored.}
}
\examples{
# Estimate first order Markov model based on the sequence provided
# in seq.fasta

library(mdist)
seqfile=system.file("extdata","seq.fasta", package="mdist")
bgpath=system.file("extdata", package="mdist")
##bgfile=sprintf("%s/background.bgbin",bgpath)
bgfile=paste(bgpath,"background.bgbin", collapse="", sep="/")
read.background(file=seqfile,1)
store.background(file=bgfile)
read.background(file=bgfile)
}
