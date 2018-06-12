context("Background")

test_that("readBackground", {
    # check if all the files can be loaded without
    # syntax errors
    # this file contains variable length sequences and 'N's
    seqfile = system.file("extdata", "seq1.fasta", package = "motifcounter")
    seqs = Biostrings::readDNAStringSet(seqfile)
    
    # check if it also works with a DNAString
    bg = readBackground(seqs[[1]], 1)
    
    # Check if all background models are plausible
    # all must represent probabilities
    for (m in 0:3) {
        
        bg = readBackground(seqs, m)
        validObject(bg)
        # check if background slot dimensions are correctly checked
        bg_altered = bg
        bg_altered@station = 1
        expect_error(validObject(bg_altered))
        
        bg_altered = bg
        bg_altered@trans = 1
        expect_error(validObject(bg_altered))
        
        # check background slot normalization
        bg_altered = bg
        bg_altered@trans = rep(0, length(bg@trans))
        expect_error(validObject(bg_altered))
        
        bg_altered = bg
        bg_altered@station = rep(0, length(bg@station))
        expect_error(validObject(bg_altered))
    }

    # Check wrong order
    expect_error(readBackground(seqs, -1)) # negative order
    expect_error(readBackground(seqs, "a")) # wrong order type
    
    # Check wrong seqs class
    expect_error(readBackground(seq[[1]], 2)) # not DNAStringSet object
    
    # check short sequences
    seqs = Biostrings::DNAStringSet("a")
    
    # not all nucleotids seen
    expect_error(readBackground(Biostrings::DNAStringSet("a"), 0))
    expect_error(readBackground(Biostrings::DNAStringSet(""), 0))
    expect_error(readBackground("a", 1)) # not a DNAStringSet
    expect_error(readBackground("", 2)) # not a DNAStringSet
    
    # Check if the number of counts is correct
    seq = Biostrings::DNAStringSet("acgtaagg")
    expect_equal(readBackground(seq, 0)@counts, rep(8, 4))
    expect_equal(readBackground(seq, 1)@counts,
                 c(2, 2, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 1, 2, 2))
    
    
    # Check order-0 background for another toy sequence
    seqfile = system.file("extdata", "test.fa", package = "motifcounter")
    seqs = Biostrings::readDNAStringSet(seqfile)
    bg = readBackground(seqs, 0)
    #c+g=21
    #a+t=15
    correct = c(30, 56, 56, 30)
    expect_equal(bg@counts, c(30, 56, 56, 30))
    expect_true(all(bg@station == correct/172))
    
    # Check order-1 background for another toy sequence
    bg = readBackground(seqs, 1)
    
    #the correct number of observations
    
    expect_equal(bg@counts, c(8,  9,  9,  2,  9, 14, 18,  9,  9, 18, 14,  9,  2,  9,  9,  8))
    correct = matrix(c(8,  9,  9,  2,  9, 14, 18,  9,  9, 18, 14,  9,  2,  9,  9,  8), 4, 4)
    correct = correct / apply(correct, 1, sum)
    expect_that(all(t(matrix(bg@trans, 4, 4)) == correct), is_true())
    
})


