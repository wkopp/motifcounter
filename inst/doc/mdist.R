library(mdist)

pvalue=0.001
gran=0.01
mdist.option(pvalue, gran)

pwmname="x3.pwm"

seq=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata",pwmname, package="mdist")

#load background as first order Markov model
read.background(pwmfile,1)

#load motif model from file
read.motif(pwmfile,"tab", 0.01)


#simluate score distribution
sims=simulate.scores(1000,1000)
plot(sims[[1]],sims[[2]])


#compute exact score distribution
dp=score.dist(pvalue, gran)
points(dp[[1]],dp[[2]], col="green")

num.motifhits("mdist/extdata/cluster10.bed.fasta")

#simulates number of motif hits and returns a distribution over the counts
# this function is only used for debug purposes!
simc=simulate.counts(seqlen=150,numofseqs=100,maxhits=400,nperm=1000)

#return a distribution for the number of motif hits in a random sequences
dist=comp.pois(seqlen=150, numofseqs=100, maxhits=1000,maxclumpsize=30 )

dist$dist

# this function scans a piece of DNA and counts the number of motif hits
# based on the chosen significance level for the score threshold
obs=num.motifhits(seqfile)

obs
# this function conducts a statistical test as to whether the observed number
# of hits is significant.

p.value=comp.pois.test(obs, maxhit=150)
