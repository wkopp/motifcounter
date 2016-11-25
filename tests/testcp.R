library(mdist)

alpha=0.01
gran=0.1
seqlen=20
numofseqs=100
maxhits=200
mdistOption(alpha, gran)

for ( pwmname in c("x31.tab","x32.tab")) {
    seqfile=system.file("extdata","seq.fasta", package="mdist")
    pwmfile=system.file("extdata",pwmname, package="mdist")

    readBackground(seqfile,1)
    readMotif(pwmfile)

    op=probOverlapHit()
    seqlen=rep(100,100)

    cpdist=compoundPoissonDist(seqlen, op)

    nom=numMotifHits(seqfile)

}
