context("Profiles")

test_that("score profile", {
    alpha=0.01
    gran=0.1
    motifcounterOption(alpha, gran)

    seqfile=system.file("extdata","test2.fa", package="motifcounter")
    # palindrom
    name="x3.tab"
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))

    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)

    # error because it is no DNAStringSet
    expect_error(motifHitProfile(seqs[[1]],motif,bg)) 
    expect_error(scoreSequenceProfile(seqs[[1]],motif,bg)) 
    # error because unequal sequence length
    expect_error(motifHitProfile(seqs[1:2],motif,bg))
    expect_error(scoreSequenceProfile(seqs[1:2],motif,bg))
    #profile for palindrome
    profile=motifHitProfile(seqs[1],motif,bg)

    #check length
    expect_equal(length(profile[[1]]),length(seqs[[1]])-ncol(motif)+1)
    expect_equal(length(profile[[2]]),length(seqs[[1]])-ncol(motif)+1)

    expect_equal(profile[[1]],profile[[2]])
    x=rep(0,20)
    x[13]=1  # at this position there should be a hit
    expect_equal(profile[[1]],x)
    
    # the max score position must agree with the motif hit position
    sprofile=scoreSequenceProfile(seqs[1],motif,bg)
    expect_equal(which(sprofile$fscores==max(sprofile$fscores)),
            which(profile$fhits==1))

    # aaa repeat
    name="x8.tab"
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))


    profile=motifHitProfile(seqs[2],motif,bg)
    x=rep(0,10-4+1)
    # in the reverse strand no hit
    expect_equal(profile[[2]],x)
    x[4]=1
    # in the forward strand one hit
    expect_equal(profile[[1]],x)

    sprofile=scoreSequenceProfile(seqs[2],motif,bg)
    expect_equal(which(sprofile$fscores==max(sprofile$fscores)),
            which(profile$fhits==1))
})
