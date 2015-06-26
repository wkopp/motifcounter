\name{num.motifhits}
\alias{num.motifhits}
\title{Number of motif hits in a given DNA sequence}
\description{
The function scans a given DNA sequence and counts the number of
motif hits on both, the forward and the reverse strand.
}
\usage{
num.motifhits(seqfile)
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
noc=num.motifhits(seqfile)
noc
}

