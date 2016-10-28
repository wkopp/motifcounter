\name{simulateNumHitsDist}
\alias{simulateNumHitsDist}
\title{Simulate number of PWM hit distribution}
\description{
This function simulates DNA sequences according to the background model and
subsequently counts how many motif hits occur.
}
\usage{
simulateNumHitsDist(seqlen,maxhits,nsim)
}
\arguments{
\item{seqlen}{Integer-valued vector contaning the sequence lengths
of the individual sequences.}
\item{maxhits}{Maximal number of hits.}
\item{nsim}{Number of random samples.}
}
\examples{

library(mdist)
mdistOption(0.01, 0.01)
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata","x31.pwm",package="mdist")

readBackground(seqfile,1)
readMotif(pwmfile)

seqlen=rep(150,100)
simc=simulateNumHitsDist(seqlen,maxhits=1000,nsim=1000)
}
