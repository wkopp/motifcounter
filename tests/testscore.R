
library(mdist)
alpha=0.01
gran=0.01
seqlen=100

mdist.option(alpha, gran)
pwmname="x31.pwm"
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata",pwmname, package="mdist")

pwm=read.table(pwmfile)

read.background(seqfile)
#load motif model from file
read.motif(pwmfile,"tab", 0.1)
#simluate score distribution
sims=sim.scores(seqlen,1000)
plot(sims[[1]],sims[[2]])
#compute exact score distribution
dp=score.dist()
points(dp[[1]],dp[[2]], col="green")

