
sim.scores=function(seqlen, nsim) {
  scorerange=integer(1)
  scorerange=.C("Rscorerange",as.integer(scorerange))[[1]]
  scores=numeric(scorerange); dist=numeric(scorerange)
  length(scores)
  .C("RsimulateScores", as.numeric(scores), 
    as.numeric(dist), as.integer(seqlen), 
    as.integer(nsim))
}

sim.counts=function(seqlen, numofseqs, maxhits, nsim) {
  dist=numeric(maxhits+1);
  .C("RsimulateCountDistribution", as.numeric(dist), 
    as.integer(nsim), as.integer(seqlen),as.integer(maxhits),as.integer(numofseqs))
}

