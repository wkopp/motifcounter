library(mdist)

alpha=0.1
gran=0.1
seqlen=1000
numofseqs=1
maxhits=100


seqfile=system.file("extdata","seq.fasta", package="mdist")
read.background(seqfile,1)
print.background()
seqfile=system.file("extdata","seq1.fasta", package="mdist")
read.background(seqfile,1)
print.background()
seqfile=system.file("extdata","seq2.fasta", package="mdist")
read.background(seqfile,1)
print.background()

seqfile=system.file("extdata","cpg.fa", package="mdist")
num.sequences(seqfile)
read.background(seqfile,1)
print.background()

#bgfile=system.file("extdata","background.bgbin", package="mdist")
# num.sequences(seqfile)
# len.sequences(seqfile)


