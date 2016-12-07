context("Compound Poisson approx")

test_that("compound", {
    alpha=0.01
    gran=0.1
    maxhits=200
    motifcounterOption(alpha, gran)
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    for ( pwmname in c("x31.tab","x32.tab")) {

        motiffile=system.file("extdata",pwmname, package="motifcounter")

        seqs=Biostrings::readDNAStringSet(seqfile)
        bg=readBackground(seqs,1)
        motif=t(as.matrix(read.table(motiffile)))

        op=probOverlapHit(motif, bg,singlestranded=FALSE)
        seqlen=rep(100,100)

        expect_error(compoundPoissonDist(seqlen, op, method=0))

        expect_equal(sum(compoundPoissonDist(seqlen, op,
            method="kopp")$dist),1)
        expect_equal(sum(compoundPoissonDist(seqlen,
            op,method="pape")$dist),1)
        
        # test accuracy of the compound Poisson model
        dist=compoundPoissonDist(seqlen, op, method="kopp")$dist
        expect_equal(sum(dist*seq(0,length(dist)-1)),op$alpha*2*length(seqlen)*(seqlen[1]-ncol(motif)+1))

        op=probOverlapHit(motif, bg,singlestranded=TRUE)
        expect_equal(sum(compoundPoissonDist(seqlen,
            op, method="kopp")$dist),1)
        expect_error(compoundPoissonDist(seqlen, op,method="pape"))

        dist=compoundPoissonDist(seqlen, op, method="kopp")$dist
        expect_equal(sum(dist*seq(0,length(dist)-1)),op$alpha*length(seqlen)*(seqlen[1]-ncol(motif)+1))
    }
})
