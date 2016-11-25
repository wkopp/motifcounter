library(mdist)

alpha=0.01
gran=0.1
seqlen=100
numofseqs=10
maxhits=100
mdistOption(alpha, gran)

pwmname="x3.tab"
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata",pwmname, package="mdist")

readBackground(seqfile,1)
readMotif(pwmfile, 0.01)

op=probOverlapHit()

dist=combinatorialDist(seqlen, op)

plot(dist$dist)
