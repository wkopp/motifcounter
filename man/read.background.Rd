\name{read.background}
\alias{read.background}
\title{Estimates the background model}
\description{
This function reads the a DNA sequence in fasta format and estimates
and order m Markov model.
}
\usage{
read.background(file,order=1)
}
\arguments{
\item{file}{Name of the sequence file in fasta format.}

\item{order}{Defines the order of the Markov models used for the background.
Default: order=1.}
}
\examples{
# Estimate first order Markov model based on the sequence provided
# in seq.fasta

library(mdist)
seqfile=system.file("extdata","seq.fasta", package="mdist")
read.background(file=seqfile,1)

#bgfile=system.file("extdata","background.bgbin", package="mdist")
bgpath=system.file("extdata", package="mdist")
bgfile=paste(bgpath,"background.bgbin", collapse="", sep="/")
store.background(file=bgfile)
read.background(file=bgfile)
}
