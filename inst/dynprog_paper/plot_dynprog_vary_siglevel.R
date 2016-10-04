library(mdist)
library(seqLogo)
library(ggplot2)
library(foreach)
library(doParallel)

variants=data.frame(alpha=c(rep(0.05,4),
											rep(c(rep(.01,4),rep(0.001,4), rep(0.0001,4)),2)),
										seqlen=c(rep(100,4),rep(100,12),rep(500,12)),
										numofseqs=c(rep(2,4),rep(10,4),rep(100,4),rep(1000,4),
														rep(2,4),rep(20,4),rep(200,4)),
										name=rep(c("x32","x21","e47","sp1"),2*3+1),
										type=rep(c("tab","tab","transfac","transfac"),2*3+1))

gran=0.1
maxhits=100

# You probably want to change the directory in which the output stored
#paperdir=Sys.getenv("DPPAPER")
#figdir=sprintf("%s/figures", paperdir)
prjdirs=c("/project/Tf2motif/home/dynprog",
					"/project/Tf2motif/home/thesis")

markov=0
seqfile=system.file("extdata","seq.fasta", package="mdist")
read.background(seqfile,markov)
read.background.sampling(seqfile,markov)

registerDoParallel(20)

foreach (i = 1:nrow(variants)) %dopar% {
	alpha=variants$alpha[i]
  lalpha=-log10(alpha)
  mdist::mdist.option(alpha, gran)
  numofseqs=variants$numofseqs[i]
  seqlen=rep(variants$seqlen[i],numofseqs)
  pwmname=variants$name[i]
  type=variants$type[i]

  print (sprintf("%s:alpha=%e, slen=%d, markov=%d", pwmname,alpha,seqlen,markov))

  pwmfile=system.file("extdata",paste(pwmname,".pwm",sep="",collapse=""),
   package="mdist")

  mdist::read.motif(pwmfile,type, 0.01)
  pwm=mdist::motif2matrix()

  Smat=matrix(0,ncol=maxhits,nrow=100)
  for (isim in 1:100) {
    simc=sim.counts(seqlen,maxhits,nsim=1000)
  	Smat[isim,]=simc[[1]][1:(maxhits)]
	}
	mean_s=apply(Smat,2,mean)
	q75_s=apply(Smat,2,quantile,.75)
	q25_s=apply(Smat,2,quantile,.25)


	op=overlap.prob()
	cpdist=comp.pois(seqlen, maxhits,maxclumpsize=60, op)
	dp=dynprog.count(seqlen, maxhits, op)

	dft=data.frame(hits=rep(1:maxhits -1,4),
							prob=c(mean_s,
							cpdist[[1]][1:maxhits],
							dp[[1]][1:maxhits],
							dbinom(1:maxhits -1, sum(2*(seqlen-nrow(pwm)+1)),
												 op$alpha)),
							model=factor(c(rep("emp",maxhits),
														 rep("cpk",maxhits),
											 rep("dp",maxhits),rep("bin",maxhits))),
							q25=c(q25_s, cpdist[[1]][1:maxhits],
							             dp[[1]][1:maxhits],
							            dbinom(1:maxhits -1, sum(2*(seqlen-nrow(pwm)+1)),
							     			                          op$alpha)),
							q75=c(q75_s, cpdist[[1]][1:maxhits],
							             dp[[1]][1:maxhits],
							            dbinom(1:maxhits -1, sum(2*(seqlen-nrow(pwm)+1)),
							     			                          op$alpha)))

	df=subset(dft, hits<60)
	pk=subset(df,model=="cpk")
	pe=subset(df,model=="emp")
	pp=subset(df,model=="dp")
	pb=subset(df,model=="bin")
	qsig=which(cumsum(pe$prob)>=0.95)
	dfp=subset(dft,hits %in% qsig)

	foreach (prjdir = prjdirs) %do% {

    figdir=sprintf("%s/figures", prjdir)
		postscript(file=sprintf("%s/dp_%s_slen%d_%d_m%d.eps", figdir,pwmname,
												seqlen[1],round(1/alpha),markov),
							 width=6,height=4)
		dev.id=dev.cur()

		pl=ggplot(df,aes(x=hits,y=prob, col=model, shape=model))+
			labs(y=expression("P(X)"),x="Num. of hits")+ geom_point()+
			scale_colour_manual(name="Models",
			limits=c("emp","cpk","dp","bin"), 
			values=c("black", "blue", "red", "grey"),
												labels=c(expression('P'["E"]*"(X)"),
																expression('P'["CP"]^"K"*"(X)"),
																expression('P'["DP"]*"(X)"),
																expression('P'["Bin"]*"(X)")))+
		scale_shape_manual(name="Models",
			limits=c("emp","cpk","dp","bin"), values=c(20,19,17,3),
												labels=c(expression('P'["E"]*"(X)"),
																expression('P'["CP"]^"K"*"(X)"),
																expression('P'["DP"]*"(X)"),
																expression('P'["Bin"]*"(X)")))+
		geom_errorbar(aes(ymin=q25,ymax=q75),width=.1)+ theme_bw()+
		theme(text=element_text(size=20),
					legend.position=c(0.7,.7),
					legend.text=element_text(size=18),
					legend.key.size=unit(.4,"inches"))
		print(pl)

		dev.off(which=dev.id)
		postscript(file=sprintf("%s/dp_%s_slen%d_%d_m%d_part.eps", figdir,pwmname,
														seqlen,round(1/alpha),markov),
							 width=6,height=4)
		dev.id=dev.cur()

		pl=ggplot(dfp,aes(x=hits,y=prob, col=model, shape=model))+
			labs(y=expression("P(X)"),x="Num. of hits")+ geom_point()+
			scale_colour_manual(name="Models",
				limits=c("emp","cpk","dp","bin"), 
				values=c("black", "blue", "red", "grey"),
				labels=c(expression('P'["E"]*"(X)"),
							expression('P'["CP"]^"K"*"(X)"),
							expression('P'["DP"]*"(X)"),
							expression('P'["Bin"]*"(X)")))+
				scale_shape_manual(name="Models",
				limits=c("emp","cpk","dp","bin"), values=c(20,19,17,3),
				labels=c(expression('P'["E"]*"(X)"),
							expression('P'["CP"]^"K"*"(X)"),
							expression('P'["DP"]*"(X)"),
							expression('P'["Bin"]*"(X)")))+
			geom_errorbar(aes(ymin=q25,ymax=q75),width=.1)+ theme_bw()+
			theme(text=element_text(size=20),
						legend.position=c(0.7,.7),
						legend.text=element_text(size=18),
						legend.key.size=unit(.4,"inches"))
		print(pl)

		dev.off(which=dev.id)

		r=3

    figdir=sprintf("%s/tab", prjdir)
		write(
		paste(
		markov, alpha,
		signif(sum(abs(pe$prob-pk$prob)),digits=r),
		signif(sum(abs(pe$prob-pp$prob)),digits=r),
		signif(sum(abs(pe$prob-pb$prob)),digits=r),
		signif(sum(abs(pe$prob[qsig]-pk$prob[qsig])),digits=r),
		signif(sum(abs(pe$prob[qsig]-pp$prob[qsig])),digits=r),
		signif(sum(abs(pe$prob[qsig]-pb$prob[qsig])),digits=r),sep=" & "),
		file=sprintf("%s/dp_%s_slen%d_%d_m%d.tabpart", figdir,pwmname,
								seqlen[1],round(1/alpha),markov))

    #figdir=sprintf("%s/figures", prjdir)
		write(
		paste(
		markov, alpha,seqlen[1],
		signif(sum(abs(pe$prob-pk$prob)),digits=r),
		signif(sum(abs(pe$prob-pp$prob)),digits=r),
		signif(sum(abs(pe$prob-pb$prob)),digits=r),
		#signif(sum(abs(pe$prob[qsig]-pk$prob[qsig])),digits=r),
		#signif(sum(abs(pe$prob[qsig]-pp$prob[qsig])),digits=r),
		#signif(sum(abs(pe$prob[qsig]-pb$prob[qsig])),digits=r),
		sep=" & "),
		file=sprintf("%s/dp_%s_slen%d_%d_m%d_tv.tabpart", figdir,pwmname,
								seqlen[1],round(1/alpha),markov))

		write(
		paste(
		markov, alpha,seqlen[1],
		#signif(sum(abs(pe$prob-pk$prob)),digits=r),
		#signif(sum(abs(pe$prob-pp$prob)),digits=r),
		#signif(sum(abs(pe$prob-pb$prob)),digits=r),
		signif(sum(abs(pe$prob[qsig]-pk$prob[qsig])),digits=r),
		signif(sum(abs(pe$prob[qsig]-pp$prob[qsig])),digits=r),
		signif(sum(abs(pe$prob[qsig]-pb$prob[qsig])),digits=r),sep=" & "),
		file=sprintf("%s/dp_%s_slen%d_%d_m%d_5p.tabpart", 
								figdir,pwmname,
								seqlen[1],round(1/alpha),markov))
  }
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



