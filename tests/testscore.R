
library(mdist)
alpha=0.01
gran=0.01

mdist.option(alpha, gran)
pwmname="x31.pwm"
pwmfile=system.file("extdata",pwmname, package="mdist")
read.motif(pwmfile,"tab", 0.01)
pwmname="sp1.pwm"
pwmfile=system.file("extdata",pwmname, package="mdist")
read.motif(pwmfile,"transfac", 0.01)
seqfile=system.file("extdata","seq.fasta", package="mdist")
seqfile=system.file("extdata","cpg.fa", package="mdist")

#pwm=read.table(pwmfile)
pwm=motif2matrix()
seqlen=ncol(pwm)


for (m in seq(0,4)) {
read.background(seqfile, m)
read.background.sampling(seqfile,m)
#load motif model from file
#read.motif(pwmfile,"tab", 0.1)
#simluate score distribution
sims=sim.scores(seqlen,1000000)
plot(sims[[1]],sims[[2]])
#compute exact score distribution
dp=score.dist()
points(dp[[1]],dp[[2]], col="green")
bf=score.dist.bf()
points(bf[[1]],bf[[2]], col="blue")
#points(sims[[1]],sims[[2]],col="red")
print(mean(abs(sims[[2]]-dp[[2]])))
print(mean(abs(bf[[2]]-dp[[2]])))
}

