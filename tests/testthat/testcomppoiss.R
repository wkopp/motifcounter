context("Compound Poisson approx")

test_that("compoundPoissonDist", {
    # Initial settings
    alpha=0.01
    gran=0.1
    motifcounterOptions(alpha, gran)

    # Obtain background
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)

    for ( pwmname in c("x31.tab","x32.tab")) {

        # Read motifs
        motiffile=system.file("extdata",pwmname, package="motifcounter")
        motif=t(as.matrix(read.table(motiffile)))

        # Overlapping hit probs
        op=probOverlapHit(motif, bg,singlestranded=FALSE)

        # check short sequences
        expect_equal(compoundPoissonDist(0,op)$dist[1],1) # too short sequence
        expect_equal(compoundPoissonDist(-1,op)$dist[1],1) # too short sequence
        expect_equal(compoundPoissonDist(ncol(motif)-1,op)$dist[1],1) # too short sequence

        # check wrong method 
        expect_error(compoundPoissonDist(100, op, method=0))

        # check wrong overlap
        expect_error(compoundPoissonDist(100, 1, method="kopp"))

        # check if dist sums to one
        seqlen=rep(100,100)
        expect_equal(sum(compoundPoissonDist(seqlen, op, method="kopp")$dist),1)
        expect_equal(sum(compoundPoissonDist(seqlen, op,method="pape")$dist),1)
        
        # check if the mean is correct
        dist=compoundPoissonDist(seqlen, op, method="kopp")$dist
        expect_equal(sum(dist*seq(0,length(dist)-1)),
                     op@alpha*2*length(seqlen)*(seqlen[1]-ncol(motif)+1))

        # check for single stranded scanning
        op=probOverlapHit(motif, bg,singlestranded=TRUE)

        # pape not supported
        expect_error(compoundPoissonDist(seqlen, op,method="pape"))

        # check if dist sums to one
        dist=compoundPoissonDist(seqlen, op, method="kopp")$dist
        expect_equal(sum(dist),1)

        # check if the mean is correct
        expect_equal(sum(dist*seq(0,length(dist)-1)),
                     op@alpha*length(seqlen)*(seqlen[1]-ncol(motif)+1))
    }
    
    # Check combinatorial sequence length
})

test_that("jaspar motif tests", {
    # Check if mean and normalization is correct for a large
    # set of real motifs from jaspar
    library(MotifDb)
    motifs=as.list(query(query(MotifDb,"hsapiens"),"JASPAR_CORE"))

    #  init
    motifcounterOptions()

    # Obtain background
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)

    seqlen=10000

    for (motif in motifs) {
        # normalize motifs
        motif=normalizeMotif(motif,0.001)

        # check scanning both strands
        op=probOverlapHit(motif, bg,singlestranded=FALSE)
        dist=compoundPoissonDist(seqlen, op, method="kopp")$dist
        # normalization
        expect_equal(sum(dist),1)
        # check mean
        expect_equal(sum(dist*seq(0,length(dist)-1)),
                     op@alpha*2*length(seqlen)*(seqlen[1]-ncol(motif)+1))

        # check scanning a single strand
        op=probOverlapHit(motif, bg,singlestranded=TRUE)
        dist=compoundPoissonDist(seqlen, op, method="kopp")$dist
        # normalization
        expect_equal(sum(dist),1)
        # check mean
        expect_equal(sum(dist*seq(0,length(dist)-1)),
                     op@alpha*length(seqlen)*(seqlen[1]-ncol(motif)+1))
    }

})

test_that("clump size distribution", {
    #  init
    motifcounterOptions()

    # Obtain background
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)

    motiffile = system.file("extdata", "x31.tab", package = "motifcounter")
    motif = t(as.matrix(read.table(motiffile)))
    maxclump = 20

    op = motifcounter:::probOverlapHit(motif, bg, singlestranded = TRUE)
    # single stranded not supported yet
    expect_error(motifcounter:::clumpSizeDist(maxclump, op))

    op = motifcounter:::probOverlapHit(motif, bg, singlestranded = FALSE)

    # test kopp clump size corretness
    dist = motifcounter:::clumpSizeDist(maxclump, op)$dist
    expect_equal(length(dist), maxclump)
    expect_equal(sum(dist),1)

    # test pape clump size corretness
    dist = motifcounter:::clumpSizeDist(maxclump, op, method = "pape")$dist
    expect_equal(length(dist), maxclump)
    expect_equal(sum(dist),1)

    dist = motifcounter:::simulateClumpSizeDist(motif, bg, 1000000, nsim = 1)$dist
    expect_equal(sum(dist), 1)
})
