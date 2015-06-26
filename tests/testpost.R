library(mdist)

alpha=0.01
gran=0.1
seqlen=100
numofseqs=10
maxhits=100
mdist.option(alpha, gran)

pwmname="x3.pwm"
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata",pwmname, package="mdist")

read.background(seqfile,1)
read.motif(pwmfile,"tab", 0.01)

op=overlap.prob()

dist=posterior.count(seqlen, numofseqs, maxhits, op)

dist2=posterior.count.debug(seqlen, numofseqs, maxhits, op,"poisson")
dist3=posterior.count.debug(seqlen, numofseqs, maxhits, op,"nbinom")
dist4=posterior.count.debug(seqlen, numofseqs, maxhits, op,"uniform")
dist5=posterior.count.debug(seqlen, numofseqs, maxhits, op,"truncunif")

plot(dist$dist)
points(dist2$dist,col="blue")
points(dist3$dist,col="green")
#points(dist4$dist,col="red")
points(dist5$dist,col="yellow")

p.value=comp.pois.test(num.motifhits(seqfile), op, maxhit=150)
