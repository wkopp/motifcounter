context("Background")

test_that("Background Syntax check", {
    # check if all the files can be loaded without 
    # syntax check
    for ( m in 0:3) {
        seqfile=system.file("extdata","seq.fasta", package="motifcounter")
        readBackground(seqfile,m)

        seqfile=system.file("extdata","seq1.fasta", package="motifcounter")
        readBackground(seqfile,m)

    }
})

test_that("Background correctness", {
    seqfile=system.file("extdata","test.fa", package="motifcounter")
    readBackground(seqfile,0)
    ret=fetchBackground()
    #c+g=21
    #a+t=9
    correct=c(12,21,21,12)/(66)
    expect_that(all(ret[[1]]==correct), is_true())

    readBackground(seqfile,1)
    ret=fetchBackground()
    #the correct number of observations
    correct=matrix(c(8,6,7,2,6,12,14,7,7,14,12,6,2,7,6,8),4,4)
    correct=correct/apply(correct,1,sum)
    expect_that(all(t(matrix(ret[[2]],4,4))==correct), is_true())

    # there are zero-value entries, an error occurs
    expect_error(readBackground(seqfile,2))
})
