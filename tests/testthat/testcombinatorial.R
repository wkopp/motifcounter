context("Combinatorial model")

test_that("combinatorial", {
    alpha=0.01
    gran=0.1
    seqlen=100
    numofseqs=10
    motifcounterOption(alpha, gran)

    pwmname="x3.tab"
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    motiffile=system.file("extdata",pwmname, package="motifcounter")

    readBackground(seqfile,1)
    motif=t(as.matrix(read.table(motiffile)))

    op=probOverlapHit(motif)

    dist=combinatorialDist(seqlen, op) # single sequence
    dist=combinatorialDist(rep(seqlen,numofseqs),op) # multiple sequences
    expect_error(combinatorialDist(30:100,op)) # variable length sequences
    expect_warning(combinatorialDist(30,op)) # sequence too short warning

    op=probOverlapHit(motif,singlestranded=TRUE)
    expect_error(combinatorialDist(seqlen,op)) # single strand not supported

})
