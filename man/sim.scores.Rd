\name{sim.scores}
\alias{sim.scores}
\title{Simlulate score distribution}
\description{
This function estimates the score distribution by simulation of DNA sequences
based on the background model.
}
\usage{
sim.scores(seqlen, nsim)
}
\arguments{
\item{seqlen}{Integer that defines the length of the sequence to be simulated.}

\item{nsim}{Defines the number of samples to be generated.}
}
\examples{
library(mdist)
mdist.option(0.01, 0.01)
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata","x31.pwm",package="mdist")

read.background(seqfile,1)
read.motif(pwmfile,"tab", 0.01)

sim.scores(seqlen= 1000,nsim=1000)
}
