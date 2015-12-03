library(mdist)

alpha=0.001
gran=0.01
seqlen=rep(150,100)
#numofseqs=100
maxhits=300
mdist.option(alpha, gran)

pwmname="x31.pwm"
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata",pwmname, package="mdist")

read.background(seqfile,1)
read.background.sampling(seqfile,1)
read.motif(pwmfile,"tab")

simc=sim.counts(seqlen,maxhits,nsim=1000)
#plot(simc[[1]], main="Distribution of the number of motif hits",
# xlab="# hits", ylab="Prob(hits)")
#plot(simc[[1]])

