library(mdist)
library(seqLogo)
library(ggplot2)
paperdir=Sys.getenv("CPPAPER")
paperdir="/project/Tf2motif/home/thesis"
figdir=sprintf("%s/figures", paperdir)
variants=data.frame(alpha=c(rep(.01,3),rep(0.001,3*3)),
										seqlen=c(rep(1000,3),rep(10000,3*3)),
										markov=c(rep(1,3),rep(1,3),rep(0,3),rep(2,3)),
										name=rep(c("x32","x21","e47"),4),
										type=rep(c("tab","tab","transfac"),4))
gran=0.1
maxhits=100

for (i in 1:nrow(variants)) {
	alpha=variants$alpha[i]
  lalpha=-log10(alpha)
  mdist.option(alpha, gran)
  seqlen=variants$seqlen[i]
  markov=variants$markov[i]
  pwmname=variants$name[i]

  print (sprintf("%s:alpha=%e, slen=%d, markov=%d", pwmname,alpha,seqlen,markov))

  seqfile=system.file("extdata","seq.fasta", package="mdist")
  read.background(seqfile,markov)
  read.background.sampling(seqfile,markov)

  pwmfile=system.file("extdata",paste(pwmname,".pwm",sep="",collapse=""),
  package="mdist")

  read.motif(pwmfile,variants$type[i], 0.01)
  pwm=motif2matrix()

  source(system.file("cp_paper","plotpwm.R", package="mdist"))
}

pwmname="x32"
pwmfile=system.file("extdata",paste(pwmname,".pwm",sep="",collapse=""),
 package="mdist")

read.motif(pwmfile,"tab", 0.01)
pwm=motif2matrix()

postscript(file=sprintf("%s/pal.eps", figdir))
seqLogo(pwm)
dev.off()

pwmname="x21"
pwmfile=system.file("extdata",paste(pwmname,".pwm",sep="",collapse=""),
 package="mdist")

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


