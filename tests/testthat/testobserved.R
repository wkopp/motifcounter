context("Observation")

test_that("observation", {
    alpha=0.01
    gran=0.1
    motifcounterOption(alpha, gran)

    # 1. test whether numSequences and lenSequences work
    seqfile=system.file("extdata","test2.fa", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    expect_error(numSequences(seqfile))
    nseq=numSequences(seqs)
    expect_equal(nseq,3)
    #if (nseq!=3) {
    #    stop(paste("extdata/test.fa contains 3
    #               sequences, but numSequences returned ", nseq))
    #}
    lseq=lenSequences(seqs)
    expect_error(lenSequences(seqfile))
    expect_equal(lseq,c(23,10,0))
    #if (!all(lseq==c(23,10,0))) {
    #    stop(paste("extdata/test.fa sequence lengths
    #               must be 23, 10 and 0, but lenSequences returned ", lseq))
    #}

    name="x3.tab"
    #seqfile=system.file("extdata","seq1.fasta", package="motifcounter")
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))

    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)

    nom=numMotifHits(seqs[[1]],motif,bg,singlestranded=TRUE) # using DNAString
    nom=numMotifHits(seqs[[1]],motif,bg,singlestranded=FALSE) #using DNAString
    # no DNAString
    expect_error(numMotifHits(seqfile,motif,bg,singlestranded=TRUE))

    nom=numMotifHits(seqs,motif,bg,singlestranded=TRUE) #Using DNAStringSet
    expect_equal(as.vector(nom$numofhits),c(1,0,0))

    nom=numMotifHits(seqs,motif,bg,singlestranded=FALSE) #Using DNAStringSet
    expect_equal(as.vector(nom$numofhits),c(2,0,0))

    name="x8.tab"
    seqfile=system.file("extdata","test2.fa", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))


    # in both cases, we expect 1 hit
    nom=numMotifHits(seqs,motif,bg,singlestranded=TRUE)
    expect_equal(as.vector(nom$numofhits),c(0,1,0))

    nom=numMotifHits(seqs,motif,bg,singlestranded=FALSE)
    expect_equal(as.vector(nom$numofhits),c(0,1,0))

})
