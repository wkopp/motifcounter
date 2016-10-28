\name{lenSequences}
\alias{lenSequences}
\title{Length of sequences provided in fasta file}
\description{
The function parses a fasta file and returns a vector containing
the length of each individual sequence contained in the fasta file.
Sequences containing 'N' or 'n' are skipped and set to length zero.
}
\usage{
lenSequences(seqfile)
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
lseq=lenSequences(seqfile)
lseq
}

