
#interface functions for motif representation
readBackground=function(file, order=1) {
   nseq=numSequences(file)
   lseq=lenSequences(file)
   dummy=.C("Rmakebg", as.character(file), as.integer(order),
    				 as.integer(nseq),as.integer(lseq))
}

fetchBackground=function() {
    stat=.Call("fetchStationBackground");
    trans=.Call("fetchTransBackground");
    return(list(stat,trans))
}

printBackground=function() {
  dummy=.C("RprintBackground");
}

deleteBackground=function() {
  dummy=.C("RdestroyBackground")
}

readBackgroundForSampling=function(file, order=1) {
	 nseq=numSequences(file)
   lseq=lenSequences(file)
  dummy=.C("RmakebgForSampling", as.character(file), as.integer(order),
  				 as.integer(nseq),as.integer(lseq))
}

printBackgroundForSampling=function() {
  dummy=.C("RprintBackgroundForSampling");
}

deleteBackgroundSampling=function() {
  dummy=.C("RdestroyBackgroundForSampling")
}
