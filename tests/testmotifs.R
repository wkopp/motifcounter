
library(motifcounter)

pwmfile=system.file("extdata","x1.tab", package="motifcounter")

# create a motif from tab format
readMotif(pwmfile,0)

mcorrect=matrix(c(.001,.001,.001,.997,.001,.001,.997,
                  .001,.001,.001,.997,.001,.001,.001,.001,.997),
                4,4)
if (motifLength()!=6) {
	stop(sprintf("motifLength() must equal 6 for x1.tab, but is instead %d",
				 motifLength()))
}

# return the PFM as R matrix
m=motif2matrix()
m=m[,1:4]
if (!(all(abs(m-mcorrect)<1e-6))) {
    stop("reading motif from tab format is wrong: error tolerance too high")
}

# delete the motif
deleteMotif()

# reallocate and load the same motif from an R matrix
readMotif(m)

# load a motif in transfac format

pwmfile=system.file("extdata","test.transfac", package="motifcounter")
readMotif(pwmfile)
m=motif2matrix()
vcorrect=c(.1,.1,.1,.7)+.01
vcorrect=vcorrect/sum(vcorrect)
if (!(all(abs(m[,1]-vcorrect)<1e-6))) {
    stop("reading motif from tab format is wrong: error tolerance too high")
}


