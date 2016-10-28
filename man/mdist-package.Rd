\name{mdist-package}
\alias{mdist-package}
\alias{mdist}
\docType{package}
\title{
Enrichment analysis of sequence motifs
mdist
}
\description{
The package provides functions to test for enrichment of 
motifs in terms of PWMs in DNA sequences of interest. The model is based
on the Compound Poisson Approximation. 
}
\details{
\tabular{ll}{
Package: \tab mdist\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-10-29\cr
License: \tab GPL-2\cr
}
%The method provides functions to determine
}
\author{
Wolfgang Kopp

Maintainer: Wolfgang Kopp <kopp@molgen.mpg.de>
}
%\references{
%Method not yet published.
%}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{Position frequency matrix, Motif enrichment, Motif score }
%#\seealso{
%#~~ Optional links to other man pages, e.g. ~~
%#~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%#}
\examples{
library(mdist)
# set motif scan options
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata","x31.pwm", package="mdist")
alpha=0.01
gran=0.01
seqlen=rep(100,1)
maxhits=100
mdistOption(alpha, gran)

# estimate background model from seqfile
readBackground(seqfile,1)
readBackgroundForSampling(seqfile,1)

# load motif model from pwmfile
readMotif(pwmfile,"tab", 0.01)

# compute score distribution by dynamic programming
dp=scoreDist()
plot(dp[[1]],dp[[2]], col="green")

# simulate the number of motif hits on random sequences
simc=simNumHitsDist(seqlen,maxhits,nsim=1000)

# Compute overlap probabilities
op=probOverlapHits()

# compute a compound poisson approximation of the motif hits distribution
cpdist=compoundPoissonDist(seqlen, maxhits,maxclumpsize=30, op)

pdist=combinatorialDist(seqlen, maxhits, op)

# plot the probability distributions
plot(simc[[1]], main="Distribution of the number of motif hits",
 xlab="# hits", ylab="Prob(hits)")
points(cpdist[[1]],col="green")
points(pdist[[1]],col="blue")
legend(70,.23,c("simulated","Comp. Poisson","Combinatorial Dist."), 
 fill=c("black","green", "blue"))

# Determines whether the motif is significantly enriched in seqfile and
# returns the corresponding pvalue.
#pval_cp=comp.pois.test(num.motifhits(seqfile), op, maxhits=150)
#pval_p=dynprog.count.test(num.motifhits(seqfile), op, maxhits=150)


}
