library(motifcounter)

for ( m in 0:3) {
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    readBackground(seqfile,m)

    printBackground()
    seqfile=system.file("extdata","seq1.fasta", package="motifcounter")
    readBackground(seqfile,m)
    printBackground()

}

library(motifcounter)
seqfile=system.file("extdata","test.fa", package="motifcounter")
readBackground(seqfile,0)
printBackground()
ret=fetchBackground()
#c+g=21
#a+t=9
correct=c(12,21,21,12)/(66)
if (all(ret[[1]]==correct)==FALSE) {
    stop("Background model for order zero is incorrect")
}

readBackground(seqfile,1)
printBackground()
ret=fetchBackground()
#the correct number of observations
correct=matrix(c(8,6,7,2,6,12,14,7,7,14,12,6,2,7,6,8),4,4)
correct=correct/apply(correct,1,sum)
if (all(t(matrix(ret[[2]],4,4))==correct)==FALSE) {
    stop("Background model for order one is incorrect")
}
