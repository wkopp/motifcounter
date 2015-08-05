library(devtools)
library(mdist)
library(seqLogo)
library(ggplot2)
install("mdist", args="--no-lock")
figdir="/project/Tf2motif/home/comppois_man/figures"
#variants=data.frame(alpha=c(rep(.01,3),rep(0.001,3*3)),
#										seqlen=c(rep(1000,3),rep(10000,3*3)),
#										markov=c(rep(1,3),rep(1,3),rep(0,3),rep(2,3)),
#										name=rep(c("x32","x21","e47"),4),
#										type=rep(c("tab","tab","transfac"),4))
alpha=0.01
lalpha=-log10(alpha)
gran=0.1
mdist.option(alpha, gran)
seqlen=10000
numofseqs=1
maxhits=500
markov=0
markovsampling=1

addlegend=FALSE

pwmname="sp1"

print (sprintf("%s:alpha=%e, slen=%d, markov=%d", pwmname,alpha,seqlen,markov))
	#alpha=variants$alpha[i]
# reproduce figure 2a,c,e use seqlen=1000, alpha=0.01, markov=1
# reproduce figure 2b,d,f use seqlen=10000, alpha=0.001, markov=1
# reproduce figure 3a,c,e use seqlen=10000, alpha=0.001, markov=0
# reproduce figure 3b,d,f use seqlen=10000, alpha=0.001, markov=2

seqfile=system.file("extdata","cpg.fa", package="mdist")
#seqfile=system.file("extdata","test.fa", package="mdist")
#seqfile=system.file("extdata","seq.fasta", package="mdist")
read.background(seqfile,markov)
print.background()
read.background.sampling(seqfile,markovsampling)
print.background.sampling()

#pwmname="x32"
pwmfile=system.file("extdata",paste(pwmname,".pwm",sep="",collapse=""),
 package="mdist")
#seqfile="dhsclusters/cluster10.bed.fasta"

#load motif model from file
read.motif(pwmfile,"transfac", 0.01)
pwm=motif2matrix()

#score.dist()
simres= sim.scores(1000,10000)

result=score.dist()

postscript(file=sprintf("%s/scoredist_sp1_d%d.eps",figdir, markovsampling))
plot(result[[1]],result[[2]], xlab="scores", 
		 ylab="Probability(Score)", col="blue", type="l")
lines(simres[[1]],simres[[2]],col="green")
dev.off()

# 434 was determined as the 0.01 significance threshold
# for order 3
# r=434:529
sum(abs(result[[2]][r]-simres[[2]][r]))/sum(abs(result[[2]][r]))

#sample counts
sc=sim.counts(seqlen, numofseqs, maxhits, 10000)
#compute scores
ov=overlap.prob()
cpc=comp.pois(seqlen, numofseqs, maxhits, 20, ov)

postscript(file=sprintf("%s/cntdist_sp1_d%d.eps",figdir,markovsampling))
plot(sc[[1]], type="l", xlab="# motif hits", ylab="Probability(# motif hits)")
lines(cpc$dist,col="green")
dev.off()

postscript(file=sprintf("%s/sp1.eps", figdir))
seqLogo(pwm)
dev.off()

#postscript(file=sprintf("%s/pal.eps", figdir))
#seqLogo(pwm)
#dev.off()

#source(system.file("cp_paper","plotpwm.R", package="mdist"))

