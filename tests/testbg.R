library(mdist)

alpha=0.1
gran=0.1
seqlen=1000
numofseqs=1
maxhits=100


seqfile=system.file("extdata","seq.fasta", package="mdist")

bgfile=system.file("extdata","background.bgbin", package="mdist")
read.background(seqfile,1)
store.background(bgfile)

read.background(bgfile)
