
library(mdist)
alpha=0.01
gran=0.01
seqlen=100

mdist.option(alpha, gran)
pwmname="x31.pwm"
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata",pwmname, package="mdist")

pwm=read.table(pwmfile)

read.background(seqfile, 1)
read.background.sampling(seqfile,1)
#load motif model from file
read.motif(pwmfile,"tab", 0.1)
#simluate score distribution
sims=sim.scores(seqlen,10000)
plot(sims[[1]],sims[[2]])
#compute exact score distribution
dp=score.dist()
points(dp[[1]],dp[[2]], col="green")

read.background(seqfile, 0)
read.background.sampling(seqfile,0)
#load motif model from file
read.motif(pwmfile,"tab", 0.1)
#simluate score distribution
sims=sim.scores(seqlen,10000)
plot(sims[[1]],sims[[2]])
#compute exact score distribution
dp=score.dist()
points(dp[[1]],dp[[2]], col="green")

read.background(seqfile, 2)
read.background.sampling(seqfile,2)
#load motif model from file
read.motif(pwmfile,"tab", 0.1)
#simluate score distribution
sims=sim.scores(seqlen,10000)
plot(sims[[1]],sims[[2]])
#compute exact score distribution
dp=score.dist()
points(dp[[1]],dp[[2]], col="green")

