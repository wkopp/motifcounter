
#interface functions for motif representation
readBackground=function(file, order=1) {
   nseq=numSequences(file)
   lseq=lenSequences(file)
   dummy=.C("mdist_makebg", as.character(file), as.integer(order),
    				 as.integer(nseq),as.integer(lseq))
}

fetchBackground=function() {
    stat=.Call("mdist_fetchStationBackground");
    trans=.Call("mdist_fetchTransBackground");
    return(list(stat,trans))
}

printBackground=function() {
  dummy=.C("mdist_printBackground");
}

deleteBackground=function() {
  dummy=.C("mdist_deleteBackground")
}

readBackgroundForSampling=function(file, order=1) {
	 nseq=numSequences(file)
   lseq=lenSequences(file)
  dummy=.C("mdist_makebgForSampling", as.character(file), as.integer(order),
  				 as.integer(nseq),as.integer(lseq))
}

printBackgroundForSampling=function() {
  dummy=.C("mdist_printBackgroundForSampling");
}

deleteBackgroundForSampling=function() {
  dummy=.C("mdist_deleteBackgroundForSampling")
}
