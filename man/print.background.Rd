\name{print.background}
\alias{print.background}
\title{Prints the current background model}
\description{
This function prints the currently loaded background model.
}
\usage{
print.background()
}
\examples{
# Estimate first order Markov model based on the sequence provided
# in seq.fasta

library(mdist)
seqfile=system.file("extdata","seq.fasta", package="mdist")
read.background(file=seqfile,1)
print.background()

}
