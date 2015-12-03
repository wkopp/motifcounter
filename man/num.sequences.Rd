\name{num.sequences}
\alias{num.sequences}
\title{Number of sequence in a provided fasta file}
\description{
The function parses a fasta file and returns the number of individual
sequence in it.
}
\usage{
num.sequences(seqfile)
}
\arguments{
\item{seqfile}{Name of the sequence file in fasta format.}

}
\examples{
# Estimate first order Markov model based on the sequence provided
# in seq.fasta

library(mdist)
seqfile=system.file("extdata","seq.fasta", package="mdist")

# scan the given sequence for the motif occurances
nseq=num.sequences(seqfile)
nseq
}

