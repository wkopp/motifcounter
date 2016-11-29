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

        readBackgroundForSampling(seqfile,1)
        seqs=Biostrings::readDNAStringSet(seqfile)
        readBackground(seqfile,1)
        motif=t(as.matrix(read.table(motiffile)))

        op=probOverlapHit(motif, singlestranded=FALSE)
        seqlen=rep(100,100)

        expect_error(compoundPoissonDist(seqlen, op, method=0))
        
        cpdist=compoundPoissonDist(seqlen, op, method="kopp")
        cpdist=compoundPoissonDist(seqlen, op,method="pape")

        op=probOverlapHit(motif, singlestranded=TRUE)
        cpdist=compoundPoissonDist(seqlen, op, method="kopp")
        expect_error(compoundPoissonDist(seqlen, op,method="pape"))


    }
})
