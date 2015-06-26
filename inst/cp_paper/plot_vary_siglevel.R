library(mdist)
library(seqLogo)
library(ggplot2)
figdir="/project/Tf2motif/home/comppois_man/figures"
variants=data.frame(alpha=c(rep(.01,3),rep(0.001,3*3)),
										seqlen=c(rep(1000,3),rep(10000,3*3)),
										markov=c(rep(1,3),rep(1,3),rep(0,3),rep(2,3)),
										name=rep(c("x32","x21","e47"),4),
										type=rep(c("tab","tab","transfac"),4))
#alpha=0.001
#lalpha=-log10(alpha)
gran=0.1
#mdist.option(alpha, gran)
#seqlen=10000
numofseqs=1
maxhits=100
#markov=0

addlegend=FALSE

for (i in 1:nrow(variants)) {
	alpha=variants$alpha[i]
  lalpha=-log10(alpha)
  mdist.option(alpha, gran)
  seqlen=variants$seqlen[i]
  markov=variants$markov[i]
  pwmname=variants$name[i]

  print (sprintf("%s:alpha=%e, slen=%d, markov=%d", pwmname,alpha,seqlen,markov))
	#alpha=variants$alpha[i]
# reproduce figure 2a,c,e use seqlen=1000, alpha=0.01, markov=1
# reproduce figure 2b,d,f use seqlen=10000, alpha=0.001, markov=1
# reproduce figure 3a,c,e use seqlen=10000, alpha=0.001, markov=0
# reproduce figure 3b,d,f use seqlen=10000, alpha=0.001, markov=2

seqfile=system.file("extdata","seq.fasta", package="mdist")
read.background(seqfile,markov)

#pwmname="x32"
pwmfile=system.file("extdata",paste(pwmname,".pwm",sep="",collapse=""),
 package="mdist")
#seqfile="dhsclusters/cluster10.bed.fasta"

#load motif model from file
read.motif(pwmfile,variants$type[i], 0.01)
pwm=motif2matrix()

#postscript(file=sprintf("%s/pal.eps", figdir))
#seqLogo(pwm)
#dev.off()

source(system.file("cp_paper","plotpwm.R", package="mdist"))
}

pwmname="x32"
pwmfile=system.file("extdata",paste(pwmname,".pwm",sep="",collapse=""),
 package="mdist")
#seqfile="dhsclusters/cluster10.bed.fasta"

#load motif model from file
read.motif(pwmfile,"tab", 0.01)
pwm=motif2matrix()

postscript(file=sprintf("%s/pal.eps", figdir))
seqLogo(pwm)
dev.off()

pwmname="x21"
pwmfile=system.file("extdata",paste(pwmname,".pwm",sep="",collapse=""),
 package="mdist")
#seqfile="dhsclusters/cluster10.bed.fasta"

#load motif model from file
read.motif(pwmfile,"tab", 0.01)
pwm=motif2matrix()

postscript(file=sprintf("%s/rep.eps",figdir))
seqLogo(pwm)
dev.off()

#source(system.file("cp_paper","plotpwm.R", package="mdist"))

pwmname="e47"
pwmfile=system.file("extdata",paste(pwmname,".pwm",sep="",collapse=""),
 package="mdist")
#seqfile="dhsclusters/cluster10.bed.fasta"

#load motif model from file
read.motif(pwmfile,"transfac", 0.01)
pwm=motif2matrix()

postscript(file=sprintf("%s/e47.eps",figdir))
seqLogo(pwm)
dev.off()

#source(system.file("cp_paper","plotpwm.R", package="mdist"))


