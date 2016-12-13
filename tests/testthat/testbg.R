context("Background")

test_that("Background Syntax check", {
    # check if all the files can be loaded without
    # syntax errors
    # this file contains variable length sequences and 'N's
    seqfile=system.file("extdata","seq1.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    
    expect_error(backgroundValid(seqs)) # not a Background object
    
    # Check if all background models are plausible
    # all must represent probabilities
    for ( m in 0:3) {
        bg=readBackground(seqs,m)
        backgroundValid(bg)
        # check if background sums to one
        expect_equal(sum(bg$station),1)
        # check if trans sums to one (for each prefix)
        expect_equal(sum(bg$trans),4^m)
    }
})

test_that("Background correctness", {
    #some simple checks
    seq=Biostrings::DNAStringSet("acgtaagg")
    expect_equal(readBackground(seq,0)$counts,rep(8,4))
    expect_equal(readBackground(seq,1)$counts,
        c(2,2,1,2,2,2,2,1,1,2,2,2,2,1,2,2))
    expect_error(readBackground(seq,2)) #zero entries
    expect_error(readBackground(seq[[1]],2)) # not DNAStringSet object
    expect_error(readBackground(Biostrings:DNAStringSet("a"),2)) #short sequence
    expect_error(readBackground(seq,-1)) # negative order

    seqfile=system.file("extdata","test.fa", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,0)
    #c+g=21
    #a+t=9
    correct=c(12,21,21,12)/(66)
    expect_equal(bg$counts,c(12,21,21,12)*2)
    expect_true(all(bg$station==correct))

    bg=readBackground(seqs,1)

    #the correct number of observations

    expect_equal(bg$counts,c(8,6,7,2,6,12,14,7,7,14,12,6,2,7,6,8))
    correct=matrix(c(8,6,7,2,6,12,14,7,7,14,12,6,2,7,6,8),4,4)
    correct=correct/apply(correct,1,sum)
    expect_that(all(t(matrix(bg$trans,4,4))==correct), is_true())

})

test_that("sequence length checks", {
  seqs=Biostrings::DNAStringSet("a")
  
  expect_error(readBackground(seqs,-1)) # negative order
  
  # Check background
  expect_error(readBackground(seqs,0)) # not all nucleotids seen
  expect_error(readBackground(Biostrings::DNAStringSet(""),0)) # too short
  expect_error(readBackground(seqs,1)) # too short
  expect_error(readBackground(seqs,2)) # too short
})