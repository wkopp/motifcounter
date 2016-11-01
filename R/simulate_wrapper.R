
simulateScoreDist=function(seqlen, nsim) {
  scorerange=integer(1)
  scorerange=.C("Rscorerange",as.integer(scorerange))[[1]]
  scores=numeric(scorerange); dist=numeric(scorerange)
  length(scores)
  return(.C("RsimulateScores", as.numeric(scores), 
    as.numeric(dist), as.integer(seqlen), 
    as.integer(nsim)))
}

simulateNumHitsDist=function(seqlen, maxhits, nsim, singlestranded=F) {
  if (length(seqlen)<=0) {
    stop("seqlen must be non-empty")
  }
  dist=numeric(maxhits+1);
  return(.C("RsimulateCountDistribution", as.numeric(dist), 
    as.integer(nsim), as.integer(length(seqlen)),
    as.integer(seqlen), as.integer(maxhits), as.integer(singlestranded)))
}

