numSequences=function(seqfile) {
  nseq=integer(1)
  res=.C("mdist_numSeqs",as.character(seqfile),as.integer(nseq))
    return (res[[2]])
}

lenSequences=function(seqfile) {
  nseq=numSequences(seqfile)
  lseq=integer(nseq)
  res=.C("mdist_lenSeqs",as.character(seqfile),as.integer(nseq),as.integer(lseq))
    return (res[[3]])
}

numMotifHits=function(seqfile, singlestranded=FALSE) {
  nseq=numSequences(seqfile)
  lseq=lenSequences(seqfile)
  noh=vector(mode="integer",length=nseq)
  res=.C("mdist_numberOfHits",as.character(seqfile),as.integer(noh),
      as.integer(nseq), as.integer(lseq),as.integer(singlestranded))
  return (list(nseq=nseq, lseq=lseq,
      numofhits=res[[2]]))
}

