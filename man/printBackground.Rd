\name{printBackground}
\alias{printBackground}
\title{Prints the current background model}
\description{
This function prints the currently loaded background model. 
}
\usage{
printBackground()
}
\examples{
# Estimate first order Markov model based on the sequence provided
# in seq.fasta

library(mdist)
seqfile=system.file("extdata","seq.fasta", package="mdist")
readBackground(file=seqfile,1)
printBackground()

}
