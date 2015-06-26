\name{score.dist}
\alias{score.dist}
\title{Score distribution}
\description{
This function computes the score distribution for the provided PWM and
background. The Score distribution is computed based on dynamic programming.
}
\usage{
score.dist()
}
%\arguments{ }
\examples{

library(mdist)
mdist.option(0.01, 0.01)
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata","x31.pwm",package="mdist")

read.background(seqfile,1)
read.motif(pwmfile,"tab", 0.01)

dp=score.dist()
}

