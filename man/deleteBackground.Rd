\name{deleteBackground}
\alias{deleteBackground}
\title{Delete  background model}
\description{
This function unloads the current background model and frees its allocated memeory.
}
\usage{
deleteBackground()
}
\examples{
# Estimate first order Markov model based on the sequence provided
# in seq.fasta

library(mdist)
seqfile=system.file("extdata","seq.fasta", package="mdist")
readBackground(file=seqfile,1)
deleteBackground()

}
