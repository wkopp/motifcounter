#figdir="/project/Tf2motif/home/bayes_man/figures"
simc=sim.counts(seqlen,numofseqs,maxhits,nsim=1000)
op=overlap.prob()
cpdist=comp.pois(seqlen, numofseqs, maxhits,maxclumpsize=60, op)
cpdist.pape=comp.pois.pape(seqlen, numofseqs, maxhits,maxclumpsize=100, op)
bayes=posterior.count(seqlen, numofseqs, maxhits, op)
#sample.mc(op$alpha,op$beta,op$beta3p,op$beta5p,seqlen-ncol(pwm)+1,numofseqs,300)
#sample.mc(9.067002e-03,op$beta,op$beta3p,op$beta5p,seqlen-ncol(pwm)+1,numofseqs,300)
#sample.mc(7.751138e-03,op$beta,op$beta3p,op$beta5p,seqlen-ncol(pwm)+1,numofseqs,300)
#sample.mc(3.052459e-03,op$beta,op$beta3p,op$beta5p,seqlen-ncol(pwm)+1,numofseqs,300)

#pdist=posterior.count(seqlen, numofseqs, maxhits, op)
#dist=comp.pois(seqlen, numofseqs, maxhits,maxclumpsize=30)
#p.value=comp.pois.test(num.motifhits(seqfile), maxhit=150)
#pape=comp.pois(seqlen, numofseqs, maxhits,maxclumpsize=30, method="pape")
#post=posteriorcount()

x=0:maxhits;
#postscript(file=sprintf("%s/%s_slen%d_10-%d_m%d.eps", figdir,pwmname,seqlen,lalpha,markov))
postscript(file=sprintf("/project/Tf2motif/home/poster/figures/pal_plot.eps", figdir,pwmname,seqlen,lalpha,markov))
plot(x,simc[[1]][x+1], main="",
 xlab="# hits", ylab="Prob(hits)")

#compound poisson distribution approximation
points(x,cpdist.pape[[1]][x+1],col="green", type="p")
#points(x,cpdist.pape[[1]][x+1],col="darkgreen", type="p")
points(x,bayes[[1]][x+1],col="red", type="p")
#points(dpois(0:maxhits, 2*(seqlen-nrow(pwm)+1)*numofseqs, alpha),col="gray")
#points(x,pdist[[1]][x+1], col="red",type="l")
points(x,dbinom(x, 2*(seqlen-ncol(pwm)+1)*numofseqs, op$alpha),col="gray",type="o")

#legend(200,.06,c("simulated","Comp. Pois.", "Bayesian", "Binomial"),
# fill=c("black","green","red","gray"), cex=1.3)
average=function(x) { return(seq(0,length(x)-1)%*% x)}
variance=function(x) { return(sum(x*(seq(0,length(x)-1)-average(x))^2))}

print(c(average(simc[[1]]),
average(bayes[[1]]),
average(cpdist[[1]]),
average(dbinom(x, 2*(seqlen-ncol(pwm)+1)*numofseqs, op$alpha)),

variance(simc[[1]]),
variance(bayes[[1]]),
variance(cpdist[[1]]),
variance(dbinom(x, 2*(seqlen-ncol(pwm)+1)*numofseqs, op$alpha))))

#plot(cumsum(simc[[1]]),cumsum(dist[[1]]))
#plot(cumsum(simc[[1]]),cumsum(pape[[1]]))
dev.off()
