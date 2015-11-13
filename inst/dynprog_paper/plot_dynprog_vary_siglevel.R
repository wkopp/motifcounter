library(mdist)
library(seqLogo)
library(ggplot2)
library(foreach)
library(doParallel)

variants=data.frame(alpha=c(rep(0.05,4),rep(c(rep(.01,4),rep(0.001,4), rep(0.0001,4)),2)),
										seqlen=c(rep(100,4),rep(100,12),rep(500,12)),
										numofseqs=c(rep(2,4),rep(10,4),rep(100,4),rep(1000,4),
																rep(2,4),rep(20,4),rep(200,4)),
										name=rep(c("x32","x21","e47","sp1"),2*3+1),
										type=rep(c("tab","tab","transfac","transfac"),2*3+1))

gran=0.1
maxhits=100

# You probably want to change the directory in which the output stored
paperdir=Sys.getenv("DPPAPER")
figdir=sprintf("%s/figures", paperdir)

  markov=1
addlegend=FALSE
  seqfile=system.file("extdata","seq.fasta", package="mdist")
  read.background(seqfile,markov)
  read.background.sampling(seqfile,markov)

#cl=makeCluster(20)
# Register cluster
registerDoParallel(28)

foreach (i = 1:nrow(variants)) %dopar% {
	alpha=variants$alpha[i]
  lalpha=-log10(alpha)
  mdist::mdist.option(alpha, gran)
  numofseqs=variants$numofseqs[i]
  seqlen=variants$seqlen[i]
  pwmname=variants$name[i]
  type=variants$type[i]

  print (sprintf("%s:alpha=%e, slen=%d, markov=%d", pwmname,alpha,seqlen,markov))

  pwmfile=system.file("extdata",paste(pwmname,".pwm",sep="",collapse=""),
   package="mdist")

  mdist::read.motif(pwmfile,type, 0.01)
  pwm=mdist::motif2matrix()

print(sprintf("%s/%s_slen%d_%d_m%d.eps", figdir,pwmname,
							                        seqlen,round(1/alpha),markov))
Smat=matrix(0,ncol=maxhits,nrow=100)
for (i in 1:100) {
  simc=sim.counts(seqlen,numofseqs,maxhits,nsim=1000)
  Smat[i,]=simc[[1]][1:(maxhits)]
}
mean_s=apply(Smat,2,mean)
median_s=apply(Smat,2,median)
q75_s=apply(Smat,2,quantile,.75)
q25_s=apply(Smat,2,quantile,.25)


op=overlap.prob()
cpdist=comp.pois(seqlen, numofseqs, maxhits,maxclumpsize=60, op)
db=dynprog.count(seqlen, numofseqs, maxhits, op)

df=data.frame(hits=1:maxhits -1,
							prob=mean_s,
							cpnew=cpdist[[1]][1:maxhits],
							db=db[[1]][1:maxhits],
							bin=dbinom(1:maxhits -1, 2*(seqlen-ncol(pwm)+1)*numofseqs,
												 op$alpha),
											q25=q25_s,q75=q75_s)
df2=df
postscript(file=sprintf("%s/%s_slen%d_%d_m%d.eps", figdir,pwmname,
												seqlen,round(1/alpha),markov))
dev.id=dev.cur()

df=subset(df2, hits<60)
Model="Sampled dist."
pl=ggplot(df,aes(x=hits,y=prob, col=Model))+
  geom_errorbar(aes(ymin=q25,ymax=q75),width=.1,color="black")
  pl=pl+geom_point()
pl<-pl+geom_point(data=df, aes(x=hits,y=db,col="Dyn. prog."), size=3)
pl<-pl+geom_point(data=df, aes(x=hits,y=cpnew,col="Comp. Poisson"))
pl<-pl+geom_line(data=df, aes(x=hits,y=bin,col="Binomial dist."))
pl<- pl+scale_colour_manual(values=c("gray", "blue", "green", "black"))
pl=pl + ylab("P(# of hits)") + 
xlab("# of hits") 
pl=pl+theme_bw()
 pl=pl+theme(panel.background=element_rect(fill='white', colour="black"), 
						legend.position="none")
#}
print(pl)

dev.off(which=dev.id)
qsig=which(cumsum(df2$prob)>=0.95)

r=3
write(
paste(paste(alpha, 
signif(sum(abs(df2$prob-df2$db)), digits=r),
signif(sum(abs(df2$prob-df2$cpnew)),digits=r),
signif(sum(abs(df2$prob-df2$bin)),digits=r), sep=" & "), "\\\\"),
file=sprintf("%s/%s_slen%d_%d_m%d_tv.tabpart", figdir,pwmname,
						 seqlen,round(1/alpha),markov))

write(
paste(paste(alpha, 
signif(sum(abs(df2$prob[qsig]-df2$db[qsig])),digits=r),
signif(sum(abs(df2$prob[qsig]-df2$cpnew[qsig])),digits=r),
signif(sum(abs(df2$prob[qsig]-df2$bin[qsig])),digits=r),sep=" & "), "\\\\"),
file=sprintf("%s/%s_slen%d_%d_m%d_5p.tabpart", figdir,pwmname,seqlen,
						 round(1/alpha),markov))

r=3
write(
paste(paste(alpha, 
signif(sum(abs(df2$prob-df2$db)),digits=r),
signif(sum(abs(df2$prob-df2$cpnew)), digits=r),
signif(sum(abs(df2$prob-df2$bin)),digits=r),
signif(sum(abs(df2$prob[qsig]-df2$db[qsig])),digits=r),
signif(sum(abs(df2$prob[qsig]-df2$cpnew[qsig])),digits=r),
signif(sum(abs(df2$prob[qsig]-df2$bin[qsig])),digits=r),sep=" & "), "\\\\"),
file=sprintf("%s/%s_slen%d_%d_m%d.tabpart", figdir,pwmname,
						  seqlen,round(1/alpha),markov))

average=function(x) { return(seq(0,length(x)-1)%*% x)}
variance=function(x) { return(sum(x*(seq(0,length(x)-1)-average(x))^2))}

r=3
write(
paste(paste(alpha, 
signif(average(df2$prob),digits=r),
signif(average(df2$db),digits=r),
signif(average(df2$cpnew),digits=r),
signif(average(df2$bin),digits=r),
signif(variance(df2$prob),digits=r),
signif(variance(df2$db),digits=r),
signif(variance(df2$cpnew),digits=r),
signif(variance(df2$bin),digits=r) ,sep=" & "), "\\\\"),
file=sprintf("%s/%s_%d_%d.mu_sig.tabpart", figdir,pwmname,round(1/alpha),markov))

}

pwmname="x32"
pwmfile=system.file("extdata",paste(pwmname,".pwm",sep="",collapse=""),
 package="mdist")

#load motif model from file
read.motif(pwmfile,"tab", 0.01)
pwm=motif2matrix()

postscript(file=sprintf("%s/pal.eps",figdir))
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

pwmname="sp1"
pwmfile=system.file("extdata",paste(pwmname,".pwm",sep="",collapse=""),
 package="mdist")
#seqfile="dhsclusters/cluster10.bed.fasta"

#load motif model from file
read.motif(pwmfile,"transfac", 0.01)
pwm=motif2matrix()

postscript(file=sprintf("%s/sp1.eps",figdir))
seqLogo(pwm)
dev.off()



