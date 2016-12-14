context("Simulation")

test_that("generateDNAString", {
    # Load background
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)

    # no background object
    expect_error(generateDNAString(1,1))

    # too short sequence
    expect_warning(generateDNAString(-1,bg))
    expect_warning(generateDNAString(bg$order-1,bg))
    expect_equal(length(generateDNAString(-1,bg)),0)
    expect_equal(length(generateDNAString(bg$order-1,bg)),0)

    # Normal sequence
    expect_equal(length(generateDNAString(10,bg)),10)
})

test_that("generateDNAStringSet", {
    # Load background
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)

    # no background object
    expect_error(generateDNAStringSet(1,1))

    # too short sequence
    expect_warning(generateDNAStringSet(-1,bg))
    expect_warning(generateDNAStringSet(bg$order-1,bg))
    expect_equal(length(generateDNAStringSet(-1,bg)[[1]]),0)
    expect_equal(length(generateDNAStringSet(bg$order-1,bg)[[1]]),0)

    # Normal sequences
    expect_equal(length(generateDNAStringSet(10,bg)),1)
    expect_equal(length(generateDNAStringSet(rep(10,10),bg)),10)
})

test_that("simulateNumHitsDist", {
    # Load sequences
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load background
    bg=readBackground(seqs,1)

    # Load motif
    motiffile=system.file("extdata","x31.tab",package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    for (ss in c(TRUE,FALSE)) {
        # normal set of sequences
        ret=simulateNumHitsDist(motif,bg,rep(10,10),nsim=100,singlestranded=ss)
        expect_equal(sum(ret$dist),1)

        # one sequence
        ret=simulateNumHitsDist(motif,bg,10,nsim=100,singlestranded=ss)
        expect_equal(sum(ret$dist),1)

        # too short sequences
        ret=simulateNumHitsDist(motif,bg,0,nsim=100,singlestranded=ss)
        expect_equal(sum(ret$dist),1)
        expect_equal(ret$dist[1],1)

        # too short sequences
        ret=simulateNumHitsDist(motif,bg,ncol(motif)-1,
                                nsim=100,singlestranded=ss)
        expect_equal(sum(ret$dist),1)
        expect_equal(ret$dist[1],1)
    }
    
    # test short sequences
    expect_error(simulateNumHitsDist(motif,bg,rep(10,10),nsim=0))
    expect_error(simulateNumHitsDist(motif,bg,rep(10,10),nsim=-1))
})

test_that("scoreDistEmpirical", {
    # Load sequences
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load background
    bg=readBackground(seqs,1)
    generateDNAString(100,bg) # after background loading everything should work

    # Load motif
    motiffile=system.file("extdata","x31.tab",package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    # test if it sums to one
    expect_equal(sum(scoreDistEmpirical(motif,bg,100,10)[[2]]),1)
    expect_equal(length(scoreDistEmpirical(motif,bg,100,10)[[2]]),
                 length(scoreDistEmpirical(motif,bg,100,10)[[1]]))

    # Check too short sequence
    expect_true(is.na(sum(scoreDistEmpirical(motif,bg,0,10)[[2]])))
    expect_true(is.na(sum(scoreDistEmpirical(motif,bg,ncol(motif)-1,10)[[2]])))

    #expect warning due to short sequence length
    expect_warning(scoreDistEmpirical(motif,bg,0,10))
    expect_warning(scoreDistEmpirical(motif,bg,ncol(motif)-1,10))
})
