\name{comp.pois.singlestranded}
\alias{comp.pois.singlestranded}
\title{Compound Poisson Approximation}
\description{
This function computes the compound Poisson approximation for the chosen
PWM and background model assuming a single stranded DNA or RNA sequence.
}
\usage{
comp.pois.singlestranded(seqlen, maxhits,maxclumpsize, overlap)
}
\arguments{
\item{seqlen}{Is an integer-valued vector where the i^th element determines
the length of i^th sequence to be scanned.}
\item{maxhits}{Integer that sets the maximal number of motif hits.}
\item{maxclumpsize}{Integer that defines the maximal size of a clump.}
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

#Compute overlapping probabilities
op=overlap.prob.singlestranded()

# Computes the distribution of the number of motif hits
seqlen=rep(150,100)
dist=comp.pois.singlestranded(seqlen, maxhits=1000,maxclumpsize=30, op)
plot(0:maxhits, dist$dist)
}

