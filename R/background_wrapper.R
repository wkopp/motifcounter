
#interface functions for motif representation
read.background=function(file, order=1) {
	 nseq=num.sequences(file)
   lseq=len.sequences(file)
    dummy=.C("Rmakebg", as.character(file), as.integer(order),
    				 as.integer(nseq),as.integer(lseq))
}

print.background=function() {
  dummy=.C("RprintBackground");
}

delete.background=function() {
  dummy=.C("RdestroyBackground")
}

read.background.sampling=function(file, order=1) {
	 nseq=num.sequences(file)
   lseq=len.sequences(file)
  dummy=.C("RmakebgForSampling", as.character(file), as.integer(order),
  				 as.integer(nseq),as.integer(lseq))
}

print.background.sampling=function() {
  dummy=.C("RprintBackgroundForSampling");
}

delete.background.sampling=function() {
  dummy=.C("RdestroyBackgroundForSampling")
}
