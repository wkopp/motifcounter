context("Compound Poisson approx")

test_that("compound", {
    alpha=0.01
    gran=0.1
    seqlen=20
    numofseqs=100
    maxhits=200
    motifcounterOption(alpha, gran)

    for ( pwmname in c("x31.tab","x32.tab")) {
        seqfile=system.file("extdata","seq.fasta", package="motifcounter")
        motiffile=system.file("extdata",pwmname, package="motifcounter")

        seqs=Biostrings::readDNAStringSet(seqfile)
        readBackground(seqfile,1)
        motif=t(as.matrix(read.table(motiffile)))

        op=probOverlapHit(motif)
        seqlen=rep(100,100)

        cpdist=compoundPoissonDist(seqlen, op)

        nom=numMotifHits(motif,seqs)
    }
})
