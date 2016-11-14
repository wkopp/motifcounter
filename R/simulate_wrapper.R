#' Simlulate score distribution
#' 
#' This function estimates the score distribution by simulation of DNA
#' sequences based on the background model and scanning with the PWM.
#' 
#' 
#' @param seqlen Integer that defines the length of the sequence to be
#' simulated.
#' @param nsim Defines the number of samples to be generated.
#' @examples
#' 
#' library(mdist)
#' # Set the the significance level and the score granularity
#' mdistOption(alpha=0.01, gran=0.1)
#' 
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' pwmfile=system.file("extdata","x31.tab",package="mdist")
#' 
#' # Load an order-1 background model
#' readBackground(seqfile,1)
#' 
#' # Load a motif from the pwmfile
#' readMotif(pwmfile)
#' 
#' # generate the simulated score distribution on sequences of length 1kb using 1000 samples 
#' simulateScoreDist(seqlen= 1000,nsim=1000)
#' 
#' @export
simulateScoreDist=function(seqlen, nsim) {
  scorerange=integer(1)
  scorerange=.C("mdist_scorerange",as.integer(scorerange),PACKAGE="mdist")[[1]]
  scores=numeric(scorerange); dist=numeric(scorerange)
  length(scores)
  return(.C("mdist_simulateScores", as.numeric(scores), 
    as.numeric(dist), as.integer(seqlen), 
    as.integer(nsim),PACKAGE="mdist"))
}



#' Simulate number of PWM hit distribution
#' 
#' This function simulates DNA sequences according to the background model and
#' subsequently counts how many motif hits occur.
#' 
#' 
#' @param seqlen Integer-valued vector contaning the sequence lengths of the
#' individual sequences.
#' @param maxhits Maximal number of hits.
#' @param nsim Number of random samples.
#' @param singlestranded Boolian flag that indicates whether a single strand or
#' both strands shall be scanned for motif hits.
#' @examples
#' 
#' 
#' library(mdist)
#' mdistOption(0.01, 0.01)
#' seqfile=system.file("extdata","seq.fasta", package="mdist")
#' pwmfile=system.file("extdata","x31.tab",package="mdist")
#' 
#' readBackground(seqfile,1)
#' readMotif(pwmfile)
#' 
#' seqlen=rep(150,100)
#' simc=simulateNumHitsDist(seqlen,maxhits=1000,nsim=1000,singlestranded=FALSE)
#' 
#' simc=simulateNumHitsDist(seqlen,maxhits=1000,nsim=1000,singlestranded=TRUE)
#' 
#' @export
simulateNumHitsDist=function(seqlen, maxhits, nsim, singlestranded=F) {
  if (length(seqlen)<=0) {
    stop("seqlen must be non-empty")
  }
  dist=numeric(maxhits+1);
  return(.C("mdist_simulateCountDistribution", as.numeric(dist), 
    as.integer(nsim), as.integer(length(seqlen)),
    as.integer(seqlen), as.integer(maxhits), as.integer(singlestranded),PACKAGE="mdist"))
}

