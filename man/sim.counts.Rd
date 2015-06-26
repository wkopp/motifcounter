\name{sim.counts}
\alias{sim.counts}
\title{Simulate number of PWM hit distribution}
\description{
This function simulates DNA sequences according to the background model and
subsequently counts how many motif hits occur.
}
\usage{
sim.counts(seqlen,numofseqs,maxhits,nsim)
}
\arguments{
\item{seqlen}{Sequence length.}
\item{numofseqs}{Number of sequences.}
\item{maxhits}{Maximal number of hits.}
\item{nsim}{Number of random samples.}
}
\examples{

library(mdist)
mdist.option(0.01, 0.01)
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata","x31.pwm",package="mdist")

read.background(seqfile,1)
read.motif(pwmfile,"tab", 0.01)

simc=sim.counts(seqlen=150,numofseqs=100,maxhits=1000,nsim=1000)
}
