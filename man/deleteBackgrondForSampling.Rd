\name{delete.background.sampling}
\alias{delete.background.sampling}
\title{Delete the current background model for sampling}
\description{
This function unloads the current background model for sampling.
}
\usage{
delete.background.sampling()
}
\examples{
# Estimate first order Markov model based on the sequence provided
# in seq.fasta

library(mdist)
seqfile=system.file("extdata","seq.fasta", package="mdist")
read.background.sampling(file=seqfile,1)
delete.background.sampling()

}
