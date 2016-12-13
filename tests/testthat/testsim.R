context("Simulation")

test_that("simulate numhits", {
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)

    motiffile=system.file("extdata","x31.tab",package="motifcounter")

    # Load the motif from the motiffile
    motif=t(as.matrix(read.table(motiffile)))

    ret=simulateNumHitsDist(motif,bg,rep(10,10),nsim=100,singlestranded=TRUE)
    expect_equal(sum(ret$dist),1)
    ret=simulateNumHitsDist(motif,bg,rep(10,10),nsim=100,singlestranded=FALSE)
    expect_equal(sum(ret$dist),1)
    
    # test short sequences
    expect_error(simulateNumHitsDist(motif,bg,rep(10,10),nsim=0))
    expect_error(simulateNumHitsDist(motif,bg,rep(ncol(motif)-1,10),nsim=1))
})

test_that("simulate DNA", {
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)
    generateDNAString(100,bg) # after background loading everything should work

    expect_error(generateDNAString(0,bg)) # positive length
    bg=readBackground(seqs,2)
    expect_error(generateDNAString(1,bg)) # too short sequence
})
