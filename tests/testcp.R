library(mdist)

alpha=0.01
gran=0.1
seqlen=20
numofseqs=100
maxhits=200
mdist.option(alpha, gran)

pwmname="x31.pwm"
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata",pwmname, package="mdist")

read.background(seqfile,1)
read.motif(pwmfile,"tab")

op=overlap.prob()

cpdist=comp.pois(seqlen, numofseqs, maxhits,maxclumpsize=30, op)

nom=num.motifhits(seqfile)

p.value=comp.pois.test(nom, op, maxhits)

cppape=comp.pois.pape(seqlen, numofseqs, maxhits,maxclumpsize=30, op)
p.value=comp.pois.test(num.motifhits(seqfile), op, maxhits)
