context("Overlapping hit probs")

test_that("overlap probs", {
    motiffile=system.file("extdata","x31.tab", package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")

    alpha=0.001
    gran=0.1
    motifcounterOption(alpha, gran)

    # estimate background model from seqfile
    readBackground(seqfile,1)

    op=probOverlapHit(motif,singlestranded=FALSE)
    expect_equal(op$alpha,0.0008473007)
    expect_equal(op$beta[1],0.0)
    expect_equal(op$beta3p[1],1.0)
    expect_equal(op$beta5p[5],3.940588e-02)
})
