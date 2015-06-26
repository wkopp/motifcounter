
Smat=matrix(0,ncol=maxhits,nrow=100)
for (i in 1:100) {
  simc=sim.counts(seqlen,numofseqs,maxhits,nsim=1000)
  Smat[i,]=simc[[1]][1:(maxhits)]
}
mean_s=apply(Smat,2,median)
q75_s=apply(Smat,2,quantile,.75)
q25_s=apply(Smat,2,quantile,.25)
#se_s=apply(Smat,2,sd)/sqrt(10)


op=overlap.prob()
cpdist=comp.pois(seqlen, numofseqs, maxhits,maxclumpsize=60, op)
cppape=comp.pois.pape(seqlen, numofseqs, maxhits,maxclumpsize=100, op)

df=data.frame(hits=1:maxhits -1,
							prob=mean_s,
							cpnew=cpdist[[1]][1:maxhits],
							cppape=cppape[[1]][1:maxhits],
							bin=dbinom(1:maxhits -1, 2*(seqlen-nrow(pwm)+1)*numofseqs,
												 op$alpha),
											q25=q25_s,q75=q75_s)
#pdist=posterior.count(seqlen, numofseqs, maxhits, op)
#dist=comp.pois(seqlen, numofseqs, maxhits,maxclumpsize=30)
#p.value=comp.pois.test(num.motifhits(seqfile), maxhit=150)
#pape=comp.pois(seqlen, numofseqs, maxhits,maxclumpsize=30, method="pape")
#post=posteriorcount()
df2=df
postscript(file=sprintf("%s/%s_10-%d_m%d.eps", figdir,pwmname,lalpha,markov))
df=subset(df2, hits<60)
Model="Sampled dist."
pl=ggplot(df,aes(x=hits,y=prob, col=Model))+
  geom_errorbar(aes(ymin=q25,ymax=q75),width=.1,color="black")
  #geom_line()+
  pl=pl+geom_point()
pl<-pl+geom_point(data=df, aes(x=hits,y=cpnew,col="Comp. Poisson (new)"), size=3)
pl<-pl+geom_point(data=df, aes(x=hits,y=cppape,col="Comp. Poisson (Pape)"))
pl<-pl+geom_line(data=df, aes(x=hits,y=bin,col="Binomial dist."))
pl<- pl+scale_colour_manual(values=c("gray", "blue", "red", "black"))
pl=pl + ylab("P(# of hits)") + 
xlab("# of hits") 
if (addlegend==TRUE) {
pl=pl+theme(legend.text=element_text(size=16), axis.text.x=element_text(size=16),
			axis.text.y=element_text(size=16), text=element_text(size=16))
} else {
pl=pl+theme(legend.position="none")
}
print(pl)
dev.off()
qsig=which(cumsum(df2$prob)>=0.95)

r=3
write(
paste(paste(alpha,seqlen, 
signif(sum(abs(df2$prob-df2$cpnew)),digits=r),
signif(sum(abs(df2$prob-df2$cppape)), digits=r),
signif(sum(abs(df2$prob-df2$bin)),digits=r),
signif(sum(abs(df2$prob[qsig]-df2$cpnew[qsig])),digits=r),
signif(sum(abs(df2$prob[qsig]-df2$cppape[qsig])),digits=r),
signif(sum(abs(df2$prob[qsig]-df2$bin[qsig])),digits=r),sep=" & "), "\\\\"),
file=sprintf("%s/%s_10-%d_m%d.tabpart", figdir,pwmname,lalpha,markov))

average=function(x) { return(seq(0,length(x)-1)%*% x)}
variance=function(x) { return(sum(x*(seq(0,length(x)-1)-average(x))^2))}

r=3
write(
paste(paste(alpha,seqlen, 
signif(average(df2$prob),digits=r),
signif(average(df2$cpnew),digits=r),
signif(average(df2$cppape),digits=r),
signif(average(df2$bin),digits=r),
signif(variance(df2$prob),digits=r),
signif(variance(df2$cpnew),digits=r),
signif(variance(df2$cppape),digits=r),
signif(variance(df2$bin),digits=r) ,sep=" & "), "\\\\"),
file=sprintf("%s/%s_10-%d_m%d.mu_sig.tabpart", figdir,pwmname,lalpha,markov))


#x=0:maxhits;
#postscript(file=sprintf("%s/%s_10-%d_m%d.eps", figdir,pwmname,lalpha,markov))
#plot(x,simc[[1]][x+1], main="",
# xlab="# hits", ylab="Prob(hits)")

#compound poisson distribution approximation
#points(x,cpdist[[1]][x+1],col="green", type="p")
#points(x,cppape[[1]][x+1],col="red", type="p")
#points(dpois(0:maxhits, 2*(seqlen-nrow(pwm)+1)*numofseqs, alpha),col="gray")
#points(x,pdist[[1]][x+1], col="red",type="l")
#points(x,dbinom(x, 2*(seqlen-nrow(pwm)+1)*numofseqs, op$alpha),col="gray",type="o")

#print(alpha)

#legend(37,.1,c("simulated","Comp. Poisson (new)", "Comp. Poisson (Pape)", "Binomial"),
# fill=c("black","green","red","gray"), cex=1.7)

#plot(cumsum(simc[[1]]),cumsum(dist[[1]]))
#plot(cumsum(simc[[1]]),cumsum(pape[[1]]))
#dev.off()
