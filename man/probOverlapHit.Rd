\name{probOverlapHit}
\alias{probOverlapHit}
\title{Motif overlap probabilities}
\description{
This function computes the self-overlapping probabilites of the motif
for the loaded PWM and background based on dynamic programming.
}
\usage{
probOverlapHit(singlestranded=F)
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

# compute the overlap probabilities for scanning both DNA strands
op=overlap.prob(singlestranded=F)

# compute the overlap probabilities for scanning a single DNA strand
op=overlap.prob(singlestranded=T)
}

