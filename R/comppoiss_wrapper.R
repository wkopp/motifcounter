
probOverlapHit=function(singlestranded=FALSE) {
    mlen=motif.length()
    alpha=numeric(1)
    beta=numeric(mlen)
    beta3p=numeric(mlen)
    beta5p=numeric(mlen)
    gamma=numeric(3*mlen)
    if (singlestranded==TRUE) {
        res=.C("RoverlapSingleStranded", alpha, beta, beta3p, beta5p, gamma)
    } else {
        res=.C("Roverlap", alpha, beta, beta3p, beta5p, gamma)
    }
    return (list(alpha=res[[1]],beta=res[[2]],beta3p=res[[3]],beta5p=res[[4]],
  						 gamma=res[[5]]))
}

compoundPoissonDist=function(seqlen, overlap, singlestranded=TRUE,
  maxhits=1000, maxclumpsize=60, method="kopp") {
  dist=numeric(maxhits+1)
  if (method=="kopp") {
    res=.C("Rcompoundpoisson_useBeta", overlap$alpha,
        overlap$beta, overlap$beta3p, overlap$beta5p,
        as.numeric(dist), as.integer(length(seqlen)),
        as.integer(seqlen),
        as.integer(maxhits), as.integer(maxclumpsize),as.integer(singlestranded))
    dist=res[[5]]
  } else if (method=="pape") {
    if (singlestranded==TRUE) {
        error("The Pape implementation of the compound Poisson distribution works only with double stranded scanning")
    }
    res=.C("RcompoundpoissonPape_useGamma", overlap$gamma,
        as.numeric(dist), as.integer(length(seqlen)), as.integer(seqlen),
        as.integer(maxhits), as.integer(maxclumpsize))
    dist=res[[2]]
  } else {
      error("The method must be 'kopp' or 'pape'")
  }
  return (list(dist=dist))
}

