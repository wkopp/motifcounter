\name{dynprog.count}
\alias{dynprog.count}
\title{Count distribution according to dynamic programming}
\description{
This function computes an approxmiation of the 
distribution of the number of motif hits
by summing over all combinations of obtaining k hits
in a given sequence. The summation is done by employing
dynamic programming.
}
\usage{
dynprog.count(seqlen, numofseqs, maxhits, overlap)
}
\arguments{
\item{seqlen}{Is an integer that defines the length of DNA to be scanned.}
\item{numofseqs}{Is an integer that defines the number of DNA sequences of size seqlen
that shall be considered.}
\item{maxhits}{Integer that sets the maximal number of motif hits.}
\item{overlap}{Is a list that contains the overlap probabilities for the 
motif. It is computed by overlap.prob(). }
}

\examples{

library(mdist)
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata","x31.pwm", package="mdist")
alpha=0.001
gran=0.01
mdist.option(alpha, gran)

# estimate background model from seqfile
read.background(seqfile,1)

# load motif model from pwmfile
read.motif(pwmfile,"tab", 0.01)

# Compute overlap probabilies
op=overlap.prob()

# Computes the distribution of the number of motif hits
dist=dynprog.count(seqlen=10, numofseqs=1, maxhits=10, op)
}

