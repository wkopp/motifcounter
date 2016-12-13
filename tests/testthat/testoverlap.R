context("Overlapping hit probs")

test_that("overlap probs", {
    motiffile=system.file("extdata","x31.tab", package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    alpha=0.001
    gran=0.1
    motifcounterOptions(alpha, gran)

    # estimate background model from seqfile
    bg=readBackground(seqs,1)

    op=probOverlapHit(motif,bg,singlestranded=FALSE)
    expect_equal(op$alpha,0.0008473007,tolerance=1e-7)
    expect_equal(op$beta[1],0.0)
    expect_equal(op$beta3p[1],1.0)
    expect_equal(op$beta5p[5],3.940588e-02)

    op=probOverlapHit(motif,bg,singlestranded=TRUE)
    expect_equal(op$alpha,0.0008473007,tolerance=1e-7)
    expect_equal(op$beta[1],0.0)
    expect_equal(op$beta3p,rep(0,length(op$beta3p)))
    expect_equal(op$beta5p,rep(0,length(op$beta5p)))
})
