\name{compoundPoissonDist}
\alias{compoundPoissonDist}
\title{Compound Poisson Approximation}
\description{
This function computes the distribution of the number of motif hits
in a set of random DNA sequences according to a compound Poisson approximation.
}
\usage{
compoundPoissonDist(seqlen, overlap, maxhits=1000, 
  maxclumpsize=60, method="kopp")
}
\arguments{
\item{seqlen}{Vector of integers that 
  contains the individual sub-sequence lengths. E.g. seqlen=c(100,200)
  means that the set of sequences is made up 
  of a sequence with 100 bp and 200 bp. Note that for the 'combinatorial model'
  all sub-sequences must be of equal length.}
\item{overlap}{Overlapping hit probabilities. This argument is the result of evaluating
  the method 'probOverlapHit', which is essential for both, the compound Poisson model and the combinatorial model.}
\item{maxhits}{Integer that sets the maximal number of motif hits. By default,
maxhits is set to 1000 motif hits.
In case, the user expects to find more than 1000 hits by chance in the sequence,
this parameter needs to be increased.}
\item{maxclumpsize}{Integer that defines the maximal size of a clump. Default: maxclumpsize=60}
}

\examples{

library(mdist)
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata","x31.pwm", package="mdist")
alpha=0.001
gran=0.01
mdistOption(alpha, gran)

# estimate background model from seqfile
readBackground(seqfile,1)

# load motif model from pwmfile
readMotif(pwmfile,"tab", 0.01)

# compute the distribution for scanning a single DNA strand

#Compute overlapping probabilities
op=probOverlapHit(singlestranded=TRUE)

# Computes the distribution of the number of motif hits
seqlen=rep(150,100)
dist=compoundPoissonDist(seqlen, op,singlestranded=singlestranded)
plot(0:maxhits, dist$dist)

# compute the distribution for scanning both DNA strands

#Compute overlapping probabilities
op=probOverlapHit(singlestranded=FALSE)

# Computes the distribution of the number of motif hits
seqlen=rep(150,100)
dist=compoundPoissonDist(seqlen, op)
plot(0:maxhits, dist$dist)
}

