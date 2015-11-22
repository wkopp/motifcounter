num.sequences=function(seqfile) {
  nseq=integer(1)
  res=.C("RnumSeqs",as.character(seqfile),as.integer(nseq))
    return (res[[2]])
}

len.sequences=function(seqfile) {
  nseq=num.sequences(seqfile)
  lseq=integer(nseq)
  res=.C("RlenSeqs",as.character(seqfile),as.integer(nseq),as.integer(lseq))
    return (res[[3]])
}

num.motifhits=function(seqfile) {
  noh=integer(1)
  seqlen=integer(1)
  numofseqs=integer(1)
  res=.C("RnumberOfHits",as.character(seqfile),as.integer(noh),
    as.integer(seqlen), as.integer(numofseqs))
    return (list(seqlen=res[[3]], numofseqs=res[[4]],
      numofhits=res[[2]]))
}

num.motifhits.singlestranded=function(seqfile) {
  noh=integer(1)
  seqlen=integer(1)
  numofseqs=integer(1)
  res=.C("RnumberOfHitsSingleStranded",as.character(seqfile),as.integer(noh),
    as.integer(seqlen), as.integer(numofseqs))
    return (list(seqlen=res[[3]], numofseqs=res[[4]],
      numofhits=res[[2]]))
}

