context("Combinatorial model")

test_that("combinatorialDist", {
    # inital settings

    alpha=0.01
    gran=0.1
    motifcounterOptions(alpha, gran)

    # Obtain background
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)

    # Obtain motif
    pwmname="x3.tab"
    motiffile=system.file("extdata",pwmname, package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    # overlapping hit probs
    op=probOverlapHit(motif,bg)

    # Dist must sum to one
    # single seq
    expect_equal(sum(combinatorialDist(100, op)$dist),1)
    expect_true(combinatorialDist(100, op)$dist[1]<1)
    # multiple seqs
    expect_equal(sum(combinatorialDist(rep(100,10), op)$dist),1) 
    expect_true(combinatorialDist(rep(100,10), op)$dist[1]<1) 

    # Error with variable length seq
    expect_error(combinatorialDist(30:100,op)) # variable length sequences

    # Warning with short seqs
    expect_warning(combinatorialDist(29,op)) # sequence too short warning

    # Error with too short sequences
    expect_equal(combinatorialDist(ncol(motif)-1,op)$dist[1],1)
    expect_equal(combinatorialDist(-1,op)$dist[1],1)
    expect_equal(combinatorialDist(0,op)$dist[1],1)

    # Single stranded not supported
    op=probOverlapHit(motif,bg,singlestranded=TRUE)
    expect_error(combinatorialDist(10,op)) # single strand not supported

})
