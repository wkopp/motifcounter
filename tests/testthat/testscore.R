context("Score distribution")

test_that("scorethreshold", {
    # Set the significance level and granularity for the score computation
    motifcounterOption(alpha=0.01,gran=0.1)
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")

    # Load the order-1 background model from the DNA sequence
    readBackground(seqfile,1)

    motiffile=system.file("extdata","x31.tab",package="motifcounter")

    # Load the motif from the motiffile
    motif=t(as.matrix(read.table(motiffile)))

    sth=scoreThreshold(motif)
    expect_true(sth$alpha<=0.01)

    name="x8.tab"
    seqfile=system.file("extdata","test2.fa", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))

    readBackground(seqfile,0)

    motifcounterOption(alpha=0.001,gran=0.1)
    expect_error(scoreThreshold(motif))
    motifcounterOption(alpha=0.01,gran=0.1)
    scoreThreshold(motif)
})

test_that("scorehistogram", {
    # Set the significance level and granularity for the score computation
    motifcounterOption(alpha=0.01,gran=0.1)

    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    # Load the order-1 background model from the DNA sequence
    readBackground(seqfile,1)
    readBackgroundForSampling(seqfile,1)

    seq=generateDNAString(10)

    motiffile=system.file("extdata","x31.tab",package="motifcounter")

    # Load the motif from the motiffile
    motif=t(as.matrix(read.table(motiffile)))

    expect_error(scoreHistogram(motif,seqfile)) # not a DNAString[Set]

    # scoreHistogram same range as scoreDist
    dp=scoreDist(as.matrix(motif))
    sh=scoreHistogram(motif,seq)

    expect_equal(sh[[1]],dp[[1]])

    # expect non-zero entries
    expect_equal(sum(sh$frequency),10-ncol(motif)+1)

    seq=Biostrings::replaceLetterAt(seq,3,"N")
    sh=scoreHistogram(motif,seq)
    expect_equal(sum(sh$frequency),0)

    seqs=generateDNAStringSet(rep(30,100))

    sh=scoreHistogram(motif,seqs)
    expect_equal(sum(sh$frequency),(30-ncol(motif)+1)*100) # number of observed
})

test_that("scoresequence", {
    # Set the significance level and granularity for the score computation
    motifcounterOption(alpha=0.01,gran=0.1)

    # Load the order-1 background model from the DNA sequence
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    readBackground(seqfile,1)
    readBackgroundForSampling(seqfile,1)

    seq=generateDNAString(10)


    motiffile=system.file("extdata","x31.tab",package="motifcounter")

    # Load the motif from the motiffile
    motif=t(as.matrix(read.table(motiffile)))

    # Compute the score distribution
    scores=scoreSequence(motif,seq)

    # same number of score on both strands
    expect_equal(length(scores$fscores),length(scores$rscores))

    # Nscores= Nseq-Nmotif+1
    expect_equal(length(scores$fscores),10-ncol(motif)+1)

    seq=Biostrings::replaceLetterAt(seq,3,"N")
    sh=scoreHistogram(motif,seq)
    scores=scoreSequence(motif,seq)
    expect_equal(scores$fscores,scores$rscores) # both should be very small
                                                # and equal


})

test_that("score correctness", {
    alpha=0.01
    motifcounterOption(alpha)
    motifname="x1.tab"
    motiffile=system.file("extdata",motifname, package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    # 1. test: with motif of length one

    readBackground(seqfile, 0)
    readBackgroundForSampling(seqfile,0)
    # compute the score distribution for this case in R
    bg=fetchBackground()
    s=round((log(motif[,1])-log(bg[[1]]))*10)
    srange=seq(min(s),max(s))
    p=rep(0,length(srange))
    for (i in 1:length(srange)){
        p[i]=sum(bg[[2]][which(s==srange[i])])
    }

    dp=scoreDist(as.matrix(motif[,1]))

    expect_equal(sum(dp[[2]]),1)
    expect_equal(length(srange),length(dp[[1]]))
    expect_equal(srange,round(dp[[1]]*10))
    expect_equal(p,dp[[2]])


    # test whether the range is equally long
    # test whether the zero entries in the score
    # distribution overlap perfectly
    # test with the stationary probabilities of the background model only

    for (m in seq(0,3)) {
        readBackground(seqfile, m)
        readBackgroundForSampling(seqfile,m)
        if (m==0) {
            #simluate score distribution
            sims=scoreDistEmpirical(as.matrix(motif[,1]),1,10000)
            #compute exact score distribution
            dp=scoreDist(as.matrix(motif[,1]))
            bf=scoreDistBf(as.matrix(motif[,1]))

        } else {
            #simluate score distribution
            sims=scoreDistEmpirical(as.matrix(motif[,1:m]),m,10000)
            #compute exact score distribution
            dp=scoreDist(as.matrix(motif[,1:m]))
            bf=scoreDistBf(as.matrix(motif[,1:m]))

        }

        expect_equal(sum(dp[[2]]),1)
        expect_equal(length(bf[[1]]),length(dp[[1]]))
        expect_equal(bf[[1]],dp[[1]])
        expect_equal(length(bf[[2]]),length(dp[[2]]))
        expect_equal(bf[[2]],dp[[2]])
        expect_equal(length(sims[[1]]),length(dp[[1]]))
        expect_equal(sims[[1]],dp[[1]])
        expect_true(all(!xor(sims[[2]]==0,dp[[2]]==0)))
        # This test is incorrect
        # Due to rounding differences of scores collected on either
        # DNA strand, this condition might actually be wrong
        #if (dp[[2]][1]<=0 || dp[[2]][length(dp[[2]])]<=0) {
        #  stop(paste("The first and the last
        #             score entry must be greater than zero: ",m))
        #}
    }

    # test whether the range is equally long
    # test whether the zero entries in the score distribution overlap perfectly
    # test with the stationary and the transition probabilities
    for (m in seq(1,3)) {
        readBackground(seqfile, m)
        readBackgroundForSampling(seqfile,m)

        #simluate score distribution
        sims=scoreDistEmpirical(motif[,1:(m+1)],m+1,10000)
        #compute exact score distribution
        dp=scoreDist(motif[,1:(m+1)])
        bf=scoreDistBf(motif[,1:(m+1)])


        expect_equal(sum(dp[[2]]),1)
        expect_equal(length(bf[[1]]),length(dp[[1]]))
        expect_equal(bf[[1]],dp[[1]])
        expect_equal(length(bf[[2]]),length(dp[[2]]))
        expect_equal(bf[[2]],dp[[2]])
        expect_equal(length(sims[[1]]),length(dp[[1]]))
        expect_equal(sims[[1]],dp[[1]])
        expect_true(all(!xor(sims[[2]]==0,dp[[2]]==0)))
        # This test is incorrect
        # Due to rounding differences of scores collected on either
        # DNA strand, this condition might actually be wrong
        #if (dp[[2]][1]<=0 || dp[[2]][length(dp[[2]])]<=0) {
          #stop(paste("The first and the last
                     #score entry must be greater than zero: ",m))
        #}
    }
})
