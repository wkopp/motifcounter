context("Observation")

test_that("lenSequences", {
    expect_error(lenSequences(3)) # wrong class

    # check sequence length
    seqs=Biostrings::DNAStringSet(c("aaa","aNa","aR"))
    expect_equal(lenSequences(seqs),c(3,0,0))
})
