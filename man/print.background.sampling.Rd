\name{print.background.sampling}
\alias{print.background.sampling}
\title{Prints the current background model for sampling}
\description{
This function prints the currently loaded background model.
}
\usage{
print.background.sampling()
}
\examples{
# Estimate first order Markov model based on the sequence provided
# in seq.fasta

library(mdist)
seqfile=system.file("extdata","seq.fasta", package="mdist")
read.background.sampling(file=seqfile,1)
print.background.sampling()

}
