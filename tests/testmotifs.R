
library(mdist)

pwmfile=system.file("extdata","x1.tab", package="mdist")

# create a motif from tab format
readMotif(pwmfile,0)

mcorrect=matrix(c(.001,.001,.001,.997,.001,.001,.997,
                  .001,.001,.001,.997,.001,.001,.001,.001,.997),
                4,4)
# return the PFM as R matrix
m=motif2matrix()
if (!(all(abs(m-mcorrect)<1e-6))) {
    stop("reading motif from tab format is wrong: error tolerance too high")
}

# delete the motif
deleteMotif()

# reallocate and load the same motif from an R matrix
readMotif(m)

# load a motif in transfac format

pwmfile=system.file("extdata","sp1.transfac", package="mdist")
readMotif(pwmfile)
vcorrect=c(2,1,6,2)+.01
vcorrect=vcorrect/sum(vcorrect)
if (!(all(abs(m[,1]-vcorrect)<1e-6))) {
    stop("reading motif from tab format is wrong: error tolerance too high")
}
# load a motif in jaspar format
#pwmfile=system.file("extdata","sp1.transfac", package="mdist")
#readMotif(pwmfile)

# load a motif in meme format
#pwmfile=system.file("extdata","elk1.meme", package="mdist")
#readMotif(pwmfile)
