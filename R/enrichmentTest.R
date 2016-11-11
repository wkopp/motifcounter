
motifEnrichmentTest=function(obs,overlap,method="compoundpoisson") {
  if (method=="compoundPoisson") {
    dist=compoundPoissonDist(obs$lseq, overlap)
  } else if (method=="combinatorial") {
    stop("The combinatorial model is not yet available for the enrichment test")
    if (any(obs$lseq!=obs$lseq)) {
      stop("For the combinatorial model, all subsequences must be of equal length")
    }
  }
  if (sum(obs$numofhits)>length(dist$dist)-1) {
    p=0.0;
  } else {
    p=1-sum(dist$dist[1:obs$numofhits])
  }
  return (p)
}

