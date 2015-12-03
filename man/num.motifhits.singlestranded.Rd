\name{num.motifhits.singlestranded}
\alias{num.motifhits.singlestranded}
\title{Number of motif hits in single stranded RNA or DNA sequences}
\description{
The function scans a single stranded RNA or DNA sequence or a set of sequences
and counts the number of
motif hits on both, the forward and the reverse strand.
It returns the length of each sequence supplied in the fasta file
and the number of hits observed across all sequences.
}
\usage{
num.motifhits.singlestranded(seqfile)
}
\arguments{
\item{seqfile}{Name of the sequence file in fasta format.}

}
\examples{
# Estimate first order Markov model based on the sequence provided
# in seq.fasta

library(mdist)
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata","x31.pwm", package="mdist")
alpha=0.001
gran=0.01
mdist.option(alpha, gran)

read.background(seqfile,1)
read.motif(pwmfile,"tab",0.01)

# scan the given sequence for the motif occurances
noc=num.motifhits.singlestranded(seqfile)
noc
}

