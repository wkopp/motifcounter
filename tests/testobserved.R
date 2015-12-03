
library(mdist)

alpha=0.01
gran=0.1
maxhits=200
mdist.option(alpha, gran)

pwmname="x31.pwm"
seqfile=system.file("extdata","seq.fasta", package="mdist")
pwmfile=system.file("extdata",pwmname, package="mdist")

read.background(seqfile,1)
read.motif(pwmfile,"tab")
#print.motif(pwmfile,"tab")

nom=num.motifhits(seqfile)

pwmname="x31.pwm"
seqfile=system.file("extdata","seq1.fasta", package="mdist")
pwmfile=system.file("extdata",pwmname, package="mdist")

read.background(seqfile,1)
read.motif(pwmfile,"tab")
#print.motif(pwmfile,"tab")

nom=num.motifhits(seqfile)
