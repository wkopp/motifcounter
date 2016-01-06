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
  nseq=num.sequences(seqfile)
  lseq=len.sequences(seqfile)
  noh=vector(mode="integer",length=nseq)
  res=.C("RnumberOfHits",as.character(seqfile),as.integer(noh),
    as.integer(nseq), as.integer(lseq))
    return (list(nseq=nseq, lseq=lseq,
      numofhits=res[[2]]))
}

num.motifhits.singlestranded=function(seqfile) {
  nseq=num.sequences(seqfile)
  lseq=len.sequences(seqfile)
  noh=vector(mode="integer",length=nseq)
  res=.C("RnumberOfHitsSingleStranded",as.character(seqfile),as.integer(noh),
    as.integer(nseq), as.integer(lseq))
    return (list(nseq=nseq, lseq=lseq,
      numofhits=res[[2]]))
}

