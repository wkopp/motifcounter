context("Observation")

test_that("observation", {
    alpha=0.01
    gran=0.1
    motifcounterOption(alpha, gran)

    seqfile=system.file("extdata","test2.fa", package="motifcounter")
    name="x3.tab"
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))

    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)

    # 1. test whether the number and length of sequences is correct
    seqs=Biostrings::readDNAStringSet(seqfile)
    nom=numMotifHits(seqs,motif,bg,singlestranded=TRUE)
    expect_equal(nom$nseq,3)
    expect_equal(as.vector(nom$lseq),c(23,10,0))


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
