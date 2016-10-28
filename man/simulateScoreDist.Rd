\name{simulateScoreDist}
\alias{simulateScoreDist}
\title{Simlulate score distribution}
\description{
This function estimates the score distribution by simulation of DNA sequences
based on the background model and scanning with the PWM.
}
\usage{
simulateScoreDist(seqlen, nsim)
}
\arguments{
\item{seqlen}{Integer that defines the length of the sequence to be simulated.}

\item{nsim}{Defines the number of samples to be generated.}
}
\examples{
library(mdist)
# Set the the significance level and the score granularity
mdistOption(alpha=0.01, gran=0.1)

seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata","x31.pwm",package="mdist")

# Load an order-1 background model
readBackground(seqfile,1)

# Load a motif from the pwmfile
readMotif(pwmfile)

# generate the simulated score distribution on sequences of length 1kb using 1000 samples 
simulateScoreDist(seqlen= 1000,nsim=1000)
}
