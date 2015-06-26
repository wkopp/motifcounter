\name{comp.pois.pape}
\alias{comp.pois.pape}
\title{Pape's Compound Poisson Approximation}
\description{
This function computes the compound Poisson approximation according to (Pape
et al. 2008) for the chosen
PWM and background model.
}
\usage{
comp.pois.pape(seqlen, numofseqs, maxhits,maxclumpsize, overlap)
}
\arguments{
\item{seqlen}{Is an integer that defines the length of DNA to be scanned.}
\item{numofseqs}{Is an integer that defines how many DNA sequences of size seqlen
shall be considered.}
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
op=overlap.prob()

# Computes the distribution of the number of motif hits
dist=comp.pois.pape(seqlen=150, numofseqs=100, maxhits=1000,maxclumpsize=30, op)
}

