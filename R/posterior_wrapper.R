
dynprog.count=function(seqlen, numofseqs, maxhits, overlap) {
  dist=numeric(maxhits+1)
  #numofseq=integer(1)
  ret=.C("RPosteriorProbability", 
  overlap$alpha, overlap$beta, overlap$beta3p, overlap$beta5p,
  as.numeric(dist), as.integer(seqlen),
   as.integer(maxhits), as.integer(numofseqs),"truncunif")
  return(list(dist=ret[[5]]))
}
dynprog.count.debug=function(seqlen, numofseqs, maxhits, overlap,prior) {
  dist=numeric(maxhits+1)
  #numofseq=integer(1)
  ret=.C("RPosteriorProbability", 
  overlap$alpha, overlap$beta, overlap$beta3p, overlap$beta5p,
  as.numeric(dist), as.integer(seqlen),
   as.integer(maxhits), as.integer(numofseqs),as.character(prior))
  return(list(dist=ret[[5]]))
}
dynprog.count.test=function(obs, overlap, maxhits) {
  dist=posterior.count(obs$seqlen, obs$numofseqs, maxhits, overlap)
  if (obs$numofhits>maxhits) {
    p=0.0;
  } else {
    p=1-sum(dist$dist[1:obs$numofhits])
  }
  return (p)
}

markov.prob=function(overlap,N) {
	ret=.Call("getMarkovProb",as.integer(N),overlap$alpha,overlap$beta,overlap$beta3p,overlap$beta5p);
	return(ret)
}

sample.mc=function(alpha,beta,beta3p,beta5p,slen,nos,perm) {
	ret=.Call("sample_mc",
						as.numeric(alpha),beta,beta3p,beta5p,
						as.integer(slen),as.integer(nos),as.integer(perm));
}
