context("Score distribution")

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
            sims=simulateScoreDist(as.matrix(motif[,1]),1,100000)
            #compute exact score distribution
            dp=scoreDist(as.matrix(motif[,1]))
            bf=scoreDistBf(as.matrix(motif[,1]))

        } else {
            #simluate score distribution
            sims=simulateScoreDist(as.matrix(motif[,1:m]),m,100000)
            #compute exact score distribution
            dp=scoreDist(as.matrix(motif[,1:m]))
            bf=scoreDistBf(as.matrix(motif[,1:m]))

        }

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
        sims=simulateScoreDist(motif[,1:(m+1)],m+1,1000000)
        #compute exact score distribution
        dp=scoreDist(motif[,1:(m+1)])
        bf=scoreDistBf(motif[,1:(m+1)])

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
