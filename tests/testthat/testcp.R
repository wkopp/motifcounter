context("Compound Poisson approx")

test_that("compound", {
    alpha=0.01
    gran=0.1
    maxhits=200
    motifcounterOptions(alpha, gran)
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
        expect_equal(sum(dist*seq(0,length(dist)-1)),
                     op$alpha*2*length(seqlen)*(seqlen[1]-ncol(motif)+1))

        op=probOverlapHit(motif, bg,singlestranded=TRUE)
        expect_error(compoundPoissonDist(seqlen, op,method="pape"))

        dist=compoundPoissonDist(seqlen, op, method="kopp")$dist
        expect_equal(sum(dist),1)

        expect_equal(sum(dist*seq(0,length(dist)-1)),
                     op$alpha*length(seqlen)*(seqlen[1]-ncol(motif)+1))
    }
    
    # Check combinatorial sequence length
    expect_error(compoundPoissonDist(0,op)) # too short sequence
    expect_error(compoundPoissonDist(ncol(motif)-1,op)) # too short sequence
})

test_that("jaspar motif tests", {
    library(MotifDb)
    motifs=as.list(query(query(MotifDb,"hsapiens"),"JASPAR_CORE"))

    alpha=0.001
    gran=0.1
    seqlen=10000
    motifcounterOptions(alpha, gran)
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)

    for (motif in motifs) {
        motif=normalizeMotif(motif,0.001)

        # check scanning both strands
        op=probOverlapHit(motif, bg,singlestranded=FALSE)
        dist=compoundPoissonDist(seqlen, op, method="kopp")$dist
        expect_equal(sum(dist),1)
        expect_equal(sum(dist*seq(0,length(dist)-1)),
                     op$alpha*2*length(seqlen)*(seqlen[1]-ncol(motif)+1))

        # check scanning a single strand
        op=probOverlapHit(motif, bg,singlestranded=TRUE)
        dist=compoundPoissonDist(seqlen, op, method="kopp")$dist
        expect_equal(sum(dist),1)
        expect_equal(sum(dist*seq(0,length(dist)-1)),
                     op$alpha*length(seqlen)*(seqlen[1]-ncol(motif)+1))
    }

})
