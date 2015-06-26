num.motifhits=function(seqfile) {
  noh=integer(1)
  seqlen=integer(1)
  numofseqs=integer(1)
  res=.C("RnumberOfHits",as.character(seqfile),as.integer(noh),
    as.integer(seqlen), as.integer(numofseqs))
    return (list(seqlen=res[[3]], numofseqs=res[[4]],
      numofhits=res[[2]]))
}

