\name{delete.background}
\alias{delete.background}
\title{Delete current background model}
\description{
This function unloads the current background model.
}
\usage{
delete.background()
}
\examples{
# Estimate first order Markov model based on the sequence provided
# in seq.fasta

library(mdist)
seqfile=system.file("extdata","seq.fasta", package="mdist")
read.background(file=seqfile,1)
delete.background()

}
