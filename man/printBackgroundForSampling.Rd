\name{printBackgroundForSampling}
\alias{printBackgroundForSampling}
\title{Prints the current background model for sampling}
\description{
This function prints the currently loaded background model that is employed
during scanning the sampled DNA sequence. It is not used for the 
compound Poisson model and the combinatorial model
}
\usage{
printBackgroundForSampling()
}
\examples{
# Estimate first order Markov model based on the sequence provided
# in seq.fasta

library(mdist)
seqfile=system.file("extdata","seq.fasta", package="mdist")
readBackgroundForSampling(file=seqfile,1)
printBackgroundForSampling()

}
