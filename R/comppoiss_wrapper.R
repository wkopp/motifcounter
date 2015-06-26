
overlap.prob=function() {
  mlen=motif.length()
  alpha=numeric(1)
  beta=numeric(mlen)
  beta3p=numeric(mlen)
  beta5p=numeric(mlen)
  gamma=numeric(3*mlen)
  res=.C("Roverlap", alpha, beta, beta3p, beta5p, gamma)
  return (list(alpha=res[[1]],beta=res[[2]],beta3p=res[[3]],beta5p=res[[4]],
  						 gamma=res[[5]]))
}

comp.pois=function(seqlen, numofseqs, 
  maxhits, maxclumpsize, overlap) {
  dist=numeric(maxhits+1)
  res=.C("Rcompoundpoisson_useBeta", overlap$alpha,
    overlap$beta, overlap$beta3p, overlap$beta5p,
    as.numeric(dist), as.integer(seqlen),as.integer(numofseqs),
    as.integer(maxhits), as.integer(maxclumpsize))
  return (list(dist=res[[5]]))
}

comp.pois.pape=function(seqlen, numofseqs, 
  maxhits, maxclumpsize, overlap) {
  dist=numeric(maxhits+1)
  res=.C("RcompoundpoissonPape_useGamma", overlap$gamma,
    as.numeric(dist), as.integer(seqlen),as.integer(numofseqs),
    as.integer(maxhits), as.integer(maxclumpsize))
  return (list(dist=res[[2]]))
}

comp.pois.test=function(obs, overlap, maxhits, maxclumpsize=20) {
  dist=comp.pois(obs$seqlen, obs$numofseqs, maxhits, maxclumpsize, overlap)
  if (obs$numofhits>maxhits) {
    p=0.0;
  } else {
    p=1-sum(dist$dist[1:obs$numofhits])
  }
  return (p)
}
