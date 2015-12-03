library(mdist)

# significance level for identifying motif hits
alpha=0.01
# grandularity for binarizing the score range
gran=0.1
# length of each RNA to be scanned 
seqlen=20
# number of sequences
numofseqs=100
#maximum number of hits
maxhits=200

#some initialization
mdist.option(alpha, gran)

#load the background as order-1 Markov  model
pwmname="x31.pwm"
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata",pwmname, package="mdist")
order=1
read.background(seqfile,order)

# load the motif
read.motif(pwmfile,"tab")

# compute overlapping hit probabilities
op=overlap.prob.singlestranded()

# compute the compound Poisson distribution for the number of hits
cpdist=comp.pois.singlestranded(seqlen, 
																maxhits,maxclumpsize=30, op)

# observed number of motif hits
nom=num.motifhits.singlestranded(seqfile)

p.value=comp.pois.singlestranded.test(
					num.motifhits.singlestranded(seqfile), op, maxhits)

