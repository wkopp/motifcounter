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

#dist2=posterior.count.debug(seqlen, numofseqs, maxhits, op,"poisson")
#dist3=posterior.count.debug(seqlen, numofseqs, maxhits, op,"nbinom")
#dist4=posterior.count.debug(seqlen, numofseqs, maxhits, op,"uniform")
#dist5=posterior.count.debug(seqlen, numofseqs, maxhits, op,"truncunif")

plot(dist$dist)
#points(dist2$dist,col="blue")
#points(dist3$dist,col="green")
#points(dist4$dist,col="red")
#points(dist5$dist,col="yellow")

#p.value=comp.pois.test(num.motifhits(seqfile), op, maxhit=150)
