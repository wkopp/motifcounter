
Smat=matrix(0,ncol=maxhits,nrow=100)
for (isim in 1:100) {
  simc=sim.counts(seqlen,maxhits,nsim=1000)
  Smat[isim,]=simc[[1]][1:(maxhits)]
}
mean_s=apply(Smat,2,mean)
q75_s=apply(Smat,2,quantile,.75)
q25_s=apply(Smat,2,quantile,.25)


op=overlap.prob()
cpdist=comp.pois(seqlen,  maxhits,maxclumpsize=60, op)
cppape=comp.pois.pape(seqlen, maxhits,maxclumpsize=100, op)

dft=data.frame(hits=rep(1:maxhits -1,4),
							prob=c(mean_s,
							cpdist[[1]][1:maxhits],
							cppape[[1]][1:maxhits],
							dbinom(1:maxhits -1, 2*(seqlen-nrow(pwm)+1),
												 op$alpha)),
							model=factor(c(rep("emp",maxhits),
														 rep("cpk",maxhits),
											 rep("cpp",maxhits),rep("bin",maxhits))),
							q25=c(q25_s, cpdist[[1]][1:maxhits],
							             cppape[[1]][1:maxhits],
							            dbinom(1:maxhits -1, 2*(seqlen-nrow(pwm)+1),
							     			                          op$alpha)),
							q75=c(q75_s, cpdist[[1]][1:maxhits],
							             cppape[[1]][1:maxhits],
							            dbinom(1:maxhits -1, 2*(seqlen-nrow(pwm)+1),
							     			                          op$alpha)))

df=subset(dft, hits<60)
postscript(file=sprintf("%s/%s_10-%d_m%d.eps", figdir,pwmname,lalpha,markov),
					 width=6,height=4)
pl=ggplot(df,aes(x=hits,y=prob, col=model, shape=model))+
  labs(y=expression("P(X)"),x="Num. of hits")+ geom_point()+
  scale_colour_manual(name="Models",
  	 limits=c("emp","cpk","cpp","bin"), values=c("black", "blue", "red", "grey"),
  										labels=c(expression('P'["E"]*"(X)"),
  														 expression('P'["CP"]^"K"*"(X)"),
  														 expression('P'["CP"]^"P"*"(X)"),
  														 expression('P'["Bin"]*"(X)")))+
  scale_shape_manual(name="Models",
  	 limits=c("emp","cpk","cpp","bin"), values=c(20,19,17,3),
  										labels=c(expression('P'["E"]*"(X)"),
  														 expression('P'["CP"]^"K"*"(X)"),
  														 expression('P'["CP"]^"P"*"(X)"),
  														 expression('P'["Bin"]*"(X)")))+
  geom_errorbar(aes(ymin=q25,ymax=q75),width=.1)+ theme_bw()+
  theme(text=element_text(size=20),
  			legend.position=c(0.7,.7),
  			legend.text=element_text(size=18),
  			legend.key.size=unit(.4,"inches"))
print(pl)
dev.off()

pk=subset(df,model=="cpk")
pe=subset(df,model=="emp")
pp=subset(df,model=="cpp")
pb=subset(df,model=="bin")
qsig=which(cumsum(pe$prob)>=0.95)

postscript(file=sprintf("%s/%s_10-%d_m%d_part.eps", figdir,pwmname,lalpha,markov),
					 width=6,height=4)
dfp=subset(dft,hits %in% qsig)

plp=ggplot(dfp,aes(x=hits,y=prob, col=model, shape=model))+
  labs(y=expression("P(X)"),x="Num. of hits")+ geom_point()+
  scale_colour_manual(name="Models",
  	 limits=c("emp","cpk","cpp","bin"), values=c("black", "blue", "red", "grey"),
  										labels=c(expression('P'["E"]*"(X)"),
  														 expression('P'["CP"]^"K"*"(X)"),
  														 expression('P'["CP"]^"P"*"(X)"),
  														 expression('P'["Bin"]*"(X)")))+
  scale_shape_manual(name="Models",
  	 limits=c("emp","cpk","cpp","bin"), values=c(20,19,17,3),
  										labels=c(expression('P'["E"]*"(X)"),
  														 expression('P'["CP"]^"K"*"(X)"),
  														 expression('P'["CP"]^"P"*"(X)"),
  														 expression('P'["Bin"]*"(X)")))+
  geom_errorbar(aes(ymin=q25,ymax=q75),width=.1)+ theme_bw()+
  theme(text=element_text(size=20),
  			legend.position=c(0.7,.7),
  			legend.text=element_text(size=18),
  			legend.key.size=unit(.4,"inches"))

print(plp)
dev.off()
r=3

write(
paste(
markov, sprintf("$10^%d$",lalpha),
signif(sum(abs(pe$prob-pk$prob)),digits=r),
signif(sum(abs(pe$prob-pp$prob)),digits=r),
signif(sum(abs(pe$prob-pb$prob)),digits=r),
signif(sum(abs(pe$prob[qsig]-pk$prob[qsig])),digits=r),
signif(sum(abs(pe$prob[qsig]-pp$prob[qsig])),digits=r),
signif(sum(abs(pe$prob[qsig]-pb$prob[qsig])),digits=r),sep=" & "),
file=sprintf("%s/../tab_chapter2/%s_10-%d_m%d.tabpart", figdir,pwmname,lalpha,markov))

average=function(x) { return(seq(0,length(x)-1)%*% x)}
variance=function(x) { return(sum(x*(seq(0,length(x)-1)-average(x))^2))}

