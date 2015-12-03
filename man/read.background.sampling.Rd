\name{read.background.sampling}
\alias{read.background.sampling}
\title{Estimates the background model for sampling}
\description{
This function reads the a DNA sequence in fasta format and estimates
and order m Markov model.
}
\usage{
read.background.sampling(file,order=1)
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
file=system.file("extdata","seq.fasta", package="mdist")
read.background.sampling(file,1)

}
