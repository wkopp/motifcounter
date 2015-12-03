\name{overlap.prob.singlestranded}
\alias{overlap.prob.singlestranded}
\title{Motif overlap probabilities}
\description{
This function computes the self-overlapping probabilites of the motif
for the loaded PWM and background based on dynamic programming for single
stranded DNA or RNA sequences.
}
\usage{
overlap.prob.singlestranded()
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

# compute the overlap probabilities
op=overlap.prob.singlestranded()
}

