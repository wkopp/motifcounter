
library(motifcounter)

alpha=0.01
gran=0.1
maxhits=200
motifcounterOption(alpha, gran)

# 1. test whether numSequences and lenSequences work
seqfile=system.file("extdata","test.fa", package="motifcounter")
nseq=numSequences(seqfile)
if (nseq!=3) {
    stop(paste("extdata/test.fa contains 3 
               sequences, but numSequences returned ", nseq))
}
lseq=lenSequences(seqfile)
if (!all(lseq==c(23,10,0))) {
    stop(paste("extdata/test.fa sequence lengths 
               must be 23, 10 and 0, but lenSequences returned ", lseq))
}

pwmname="x31.tab"
seqfile=system.file("extdata","seq1.fasta", package="motifcounter")
pwmfile=system.file("extdata",pwmname, package="motifcounter")

readBackground(seqfile,1)
readMotif(pwmfile)

nom=numMotifHits(seqfile)

pwmname="x31.tab"
seqfile=system.file("extdata","seq1.fasta", package="motifcounter")
pwmfile=system.file("extdata",pwmname, package="motifcounter")

readBackground(seqfile,1)
readMotif(pwmfile)

nom=numMotifHits(seqfile)
