
motifEnrichmentTest=function(obs,overlap,singlestranded,method) {
  if (method=="compoundPoisson") {
    dist=compoundPoissonDist(obs$lseq, overlap, singlestranded=TRUE)
  } else if (method=="combinatorial") {
    if (any(obs$lseq!=obs$lseq)) {
      error("For the combinatorial model, all subsequences must be of equal length")
    }
  }
  if (sum(obs$numofhits)>maxhits) {
    p=0.0;
  } else {
    p=1-sum(dist$dist[1:obs$numofhits])
  }
  return (p)
}

