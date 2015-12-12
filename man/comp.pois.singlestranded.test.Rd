\name{comp.pois.singlestranded.test}
\alias{comp.pois.singlestranded.test}
\title{Statistical enrichment test for the PWM}
\description{
This method runs a statistical test based on the compound Poisson
approximation of the number of motif hits in order to decide the whether
a PWM is enriched in a single stranded sequence of interest.
}
\usage{
comp.pois.singlestranded.test(obs, overlap, maxhits, maxclumpsize=20)
}
\arguments{
\item{obs}{Is a list containing obs$nseq, obs$lseq and obs$numofhits
which correspond to the number of sequences to be scanned,
the length of each individual sequence and
the number of identified hits, respectively.}

\item{overlap}{Is a list that contains the overlap probabilities for the 
motif. It is computed by overlap.prob(). }
\item{maxhits}{Integer that sets the maximal number of motif hits.}
\item{maxclumpsize}{Integer that sets the maximum size of a clump.}
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

# Compute overlap probability
op=overlap.prob.singlestranded()

# Determines whether the motif is significantly enriched in seqfile and
# returns the corresponding pvalue.
p.value=comp.pois.singlestranded.test(
  num.motifhits.singlestranded(seqfile), op, maxhit=150)
}
