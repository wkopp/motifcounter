numSequences=function(seqfile) {
  nseq=integer(1)
  res=.C("RnumSeqs",as.character(seqfile),as.integer(nseq))
    return (res[[2]])
}

lenSequences=function(seqfile) {
  nseq=numSequences(seqfile)
  lseq=integer(nseq)
  res=.C("RlenSeqs",as.character(seqfile),as.integer(nseq),as.integer(lseq))
    return (res[[3]])
}

numMotifHits=function(seqfile, singlestraned=TRUE) {
  nseq=numSequences(seqfile)
  lseq=lenSequences(seqfile)
  noh=vector(mode="integer",length=nseq)
  if (singlestranded==TRUE) {
    res=.C("RnumberOfHitsSingleStranded",as.character(seqfile),as.integer(noh),
      as.integer(nseq), as.integer(lseq))
  } else {
    res=.C("RnumberOfHits",as.character(seqfile),as.integer(noh),
      as.integer(nseq), as.integer(lseq))
  }
  return (list(nseq=nseq, lseq=lseq,
      numofhits=res[[2]]))
}

