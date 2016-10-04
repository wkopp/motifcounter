library(mdist)
library(seqLogo)
library(ggplot2)
paperdir=Sys.getenv("CPPAPER")
#paperdir="/project/Tf2motif/home/thesis"
figdir=sprintf("%s/figures", paperdir)
variants=data.frame(alpha=c(rep(.01,3),rep(0.001,3*3)),
										seqlen=c(rep(1000,3),rep(10000,3*3)),
										markov=c(rep(1,3),rep(1,3),rep(0,3),rep(2,3)),
										name=rep(c("x32","x21","e47"),4),
										type=rep(c("tab","tab","transfac"),4))
gran=0.1
maxhits=100
noplot=TRUE
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

	# generate results
  source(system.file("cp_paper","compare_methods.R", package="mdist"))

pl=ggplot(df,aes(x=hits,y=prob, col=model, shape=model))+
  labs(y=expression("P(X)"),x="Num. of hits")+ geom_point()+
  scale_colour_manual(name="Models",
  	 limits=c("emp","cpk","cpp","bin"), values=c("black", "blue", "red", "grey"),
  										labels=c(expression('P'["E"]*"(X)"),
  														 expression('P'["CP"]^"N"*"(X)"),
  														 expression('P'["CP"]^"I"*"(X)"),
  														 expression('P'["Bin"]*"(X)")))+
  scale_shape_manual(name="Models",
  	 limits=c("emp","cpk","cpp","bin"), values=c(20,19,17,3),
  										labels=c(expression('P'["E"]*"(X)"),
  														 expression('P'["CP"]^"N"*"(X)"),
  														 expression('P'["CP"]^"I"*"(X)"),
  														 expression('P'["Bin"]*"(X)")))+
  geom_errorbar(aes(ymin=q25,ymax=q75),width=.1)+ theme_bw()+
  theme(text=element_text(size=20),
  			legend.position=c(0.7,.7),
  			legend.text=element_text(size=18),
  			legend.key.size=unit(.4,"inches"))
print(pl)
  if (noplot==TRUE) next;
  # produce plots
  figdir="/project/Tf2motif/home/thesis/figures"
  source(system.file("cp_paper","plotpwm.R", package="mdist"))
  paperdir=Sys.getenv("CPPAPER"); figdir=sprintf("%s/figures", paperdir)
  source(system.file("cp_paper","plotpwm.R", package="mdist"))

	#produce table row entries
  source(system.file("cp_paper","gentab_paper.R", package="mdist"))
  source(system.file("cp_paper","gentab_thesis.R", package="mdist"))
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


