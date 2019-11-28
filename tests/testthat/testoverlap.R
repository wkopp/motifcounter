context("Overlap")

test_that("probOverlaphit", {

    motifcounterOptions()

    # Obtain sequence
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Obtain background
    bg=readBackground(seqs,1)

    # Obtain motif
    motiffile=system.file("extdata","x31.tab", package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    th = scoreThreshold(motif, bg)

    op=probOverlapHit(motif,bg,singlestranded=FALSE)
    expect_equal(op@alpha,th$alpha,tolerance=1e-7)
    expect_equal(op@alpha,0.0008473007,tolerance=1e-7)
    expect_equal(op@alpha, op@gamma[1])
    expect_equal(op@beta[1],0.0)
    expect_equal(op@beta3p[1],1.0)
    expect_equal(op@beta5p[5],3.940588e-02)

    op=probOverlapHit(motif,bg,singlestranded=TRUE)
    expect_equal(op@alpha,th$alpha,tolerance=1e-7)
    expect_equal(op@alpha,0.0008473007,tolerance=1e-7)
    expect_equal(op@alpha, op@gamma[1])
    expect_equal(op@beta[1],0.0)
    expect_equal(op@beta3p,rep(0,length(op@beta3p)))
    expect_equal(op@beta5p,rep(0,length(op@beta5p)))

})
