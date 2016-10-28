\name{combinatorialDist}
\alias{combinatorialDist}
\title{Combinatrial distribution for the number of motif hits}
\description{
This function approxmiates  the
distribution of the number of motif hits
by summing over all combinations of obtaining k hits
in a given sequence. The summation is done by employing
dynamic programming.
}
\usage{
combinatorialDist(seqlen, overlap, singlestranded=TRUE, maxhits=100)
}
\arguments{
\item{seqlen}{Is an integer-valued vector that 
defines the lengths of the individual sequences.}
\item{overlap}{Is a list that contains the overlap probabilities for the 
motif. It is computed by overlap.prob(). }
\item{maxhits}{Boolian value that indicates if only a single DNA strand is scanned
for motif hits or if both DNA strands are scanned.
Integer that sets the maximal number of motif hits.}
\item{maxhits}{Integer that sets the maximal number of motif hits.}
}

\examples{

library(mdist)
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata","x31.pwm", package="mdist")
alpha=0.001
mdistOption(alpha)

# estimate background model from seqfile
readBackground(seqfile,1)

# load motif model from pwmfile
readMotif(pwmfile,"tab", 0.01)

# Compute overlap probabilies
op=probOverlapHit(singlestranded=TRUE)
seqlen=rep(100,1)
# Computes the distribution of the number of motif hits
dist=combinatorialDist(seqlen, op, singlestranded=TRUE)
}

