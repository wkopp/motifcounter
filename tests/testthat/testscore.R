context("Score distribution")

test_that("scoreDists", {
    # in this part we test scoreDist, scoreDistBf and scoreDistEmpirical

    alpha=0.01
    motifcounterOptions(alpha)

    # Load sequences
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load background
    bg=readBackground(seqs, 0)

    # Load motif
    motifname="x1.tab"
    motiffile=system.file("extdata",motifname, package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    # Check motif shorter than background order
    bg=readBackground(seqs, 2)
    expect_error(scoreDist(as.matrix(motif[,1]),bg))

    # Reinit background
    bg=readBackground(seqs, 0)

    # 1. test: with motif of length one
    # compute the score distribution for this case in R
    s=round((log(motif[,1])-log(bg@station))*10)
    srange=seq(min(s),max(s))
    p=rep(0,length(srange))
    for (i in 1:length(srange)){
        p[i]=sum(bg@station[which(s==srange[i])])
    }

    dp=scoreDist(as.matrix(motif[,1]),bg)

    # test for equality of the two variants
    expect_equal(sum(dp[[2]]),1)
    expect_equal(length(srange),length(dp[[1]]))
    expect_equal(srange,round(dp[[1]]*10))
    expect_equal(p,dp[[2]])


    # test with the stationary probabilities of the background model only
    for (m in seq(0,3)) {
        bg=readBackground(seqs, m)

        if (m==0) {
            #simluate score distribution
            sims=scoreDistEmpirical(as.matrix(motif[,1]), bg,1,1000)
            #compute exact score distribution
            dp=scoreDist(as.matrix(motif[,1]),bg)
            bf=scoreDistBf(as.matrix(motif[,1]),bg)

        } else {
            #simluate score distribution
            sims=scoreDistEmpirical(as.matrix(motif[,1:m]),bg,m,1000)
            #compute exact score distribution
            dp=scoreDist(as.matrix(motif[,1:m]),bg)
            bf=scoreDistBf(as.matrix(motif[,1:m]),bg)
        }

        expect_equal(sum(dp[[2]]),1)
        expect_equal(bf[[1]],dp[[1]])
        expect_equal(bf[[2]],dp[[2]])
        expect_equal(sims[[1]],dp[[1]])
        expect_equal(sum(sims[[2]][dp[[2]]==0]),0)
    }

    # test with the stationary and the transition probabilities
    for (m in seq(1,3)) {
        bg=readBackground(seqs, m)


        # simluate score distribution and
        # compute exact score distribution
        dp=scoreDist(motif[,1:(m+1)],bg)
        bf=scoreDistBf(motif[,1:(m+1)],bg)
        sims=scoreDistEmpirical(motif[,1:(m+1)],bg,m+1,1000)


        expect_equal(sum(dp[[2]]),1)
        expect_equal(bf[[1]],dp[[1]])
        expect_equal(bf[[2]],dp[[2]])

        expect_equal(sims[[1]],dp[[1]])
        # for impossible scores, check if sampling indeed
        # yields no observations
        expect_equal(sum(sims[[2]][dp[[2]]==0]),0)

    }
})

test_that("further scoreDist checks with jaspar motifs", {
    library(MotifDb)
    motifs=as.list(query(query(MotifDb,"hsapiens"),"JASPAR_CORE"))

    # Load sequences
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load background
    bg=readBackground(seqs, 1)

    for ( motif in motifs ) {
        motif=normalizeMotif(motif)
        if (ncol(motif)>4) {
            motif=motif[, 1:4]
        }
        dp=scoreDist(as.matrix(motif),bg)
        bf=scoreDistBf(as.matrix(motif),bg)

        expect_equal(dp$dist, bf$dist)
        expect_equal(sum(dp$dist),1)

        expect_true(dp$dist[1]>0)
        expect_true(dp$dist[length(dp[[2]])]>0)
    }
})

test_that("scoreSequence", {
    # This implicitly tests scoreStrand as well,
    # so no separate tests are available

    # Set the significance level and granularity for the score computation
    motifcounterOptions(alpha=0.01,gran=0.1)

    # Load sequences
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load background
    bg=readBackground(seqs,1)

    # Load motif
    motiffile=system.file("extdata","x31.tab",package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    # no motif, not background object, no DNAString
    expect_error(scoreSequence(seqs,motif,1))
    expect_error(scoreSequence(seqs,1,bg))
    expect_error(scoreSequence("a",motif,bg))

    # expect error with too short sequence
    scores=scoreSequence(generateDNAString(ncol(motif)-1,bg),motif,bg)
    expect_equal(length(scores[[1]]),0)
    expect_equal(length(scores[[2]]),0)
    scores=scoreSequence(Biostrings::DNAString(""),motif,bg)
    expect_equal(length(scores[[1]]),0)
    expect_equal(length(scores[[2]]),0)

    # check motif shorter than bg order
    bg=readBackground(seqs,2)
    expect_error(scoreSequence(seqs,as.matrix(motif[, 1]),bg)[[1]])

    # Reload background
    bg=readBackground(seqs,1)

    seq=generateDNAString(10,bg)

    # Compute score sequence
    scores=scoreSequence(seq,motif,bg)

    # same number of positions on both strands
    scores=scoreSequence(seq,motif,bg)
    expect_equal(length(scores$fscores),length(scores$rscores))

    # and length is Nscores= Nseq-Nmotif+1
    expect_equal(length(scores$fscores),10-ncol(motif)+1)

    # If sequence contains 'N', the respective scores should be NaN
    seq=Biostrings::replaceLetterAt(seq,3,"N")
    scores=scoreSequence(seq,motif,bg)
    expect_equal(scores$fscores,scores$rscores) # both should be equal
    expect_equal(scores$fscores, rep(NaN, 3))
    expect_equal(length(scores$fscores),3) # length=3
    expect_equal(length(scores$rscores),3) # length=3
})

test_that("scoreProfile", {
    # init
    alpha=0.01
    gran=0.1
    motifcounterOptions(alpha, gran)

    # Load seqs
    seqfile=system.file("extdata","test2.fa", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load Background
    bg=readBackground(seqs,1)

    # palindrom
    name="x3.tab"
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))

    # error because it is no DNAStringSet
    # This was updated. It used to only accept DNAStringSet,
    # but now, DNAStringSet and DNAString is accepted.
    scoreProfile(seqs[[1]],motif,bg)

    # error because unequal sequence length
    expect_error(scoreProfile(seqs[1:2],motif,bg))

    # Check too short sequence
    expect_equal(length(scoreProfile(
                Biostrings::DNAStringSet(""),motif,bg)[[1]]),0)
    expect_equal(length(scoreProfile(
                Biostrings::DNAStringSet(""),motif,bg)[[2]]),0)
    expect_equal(length(scoreProfile(
                generateDNAStringSet(ncol(motif)-1,bg),motif,bg)[[1]]),0)
    expect_equal(length(scoreProfile(
                generateDNAStringSet(ncol(motif)-1,bg),motif,bg)[[2]]),0)

    # check motif shorter than bg order
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs2=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs2,2)
    expect_error(scoreProfile(seqs,as.matrix(motif[, 1]),bg)[[1]])

    # Reload background
    bg=readBackground(seqs,1)

    # Test correct hit positions
    mh=scoreProfile(seqs[1], motif,bg)
    # because it is a palindrome
    expect_equal(mh[[1]],mh[[2]])
    expect_equal(which(mh[[1]]==max(mh[[1]])),13)

    # no hit in the second sequence
    mh=scoreProfile(seqs[2], motif,bg)
    expect_equal(mh[[1]],mh[[2]])
    expect_equal(length(mh[[1]]),length(seqs[[2]])-ncol(motif)+1)

    # NaN due to 'N'
    mh=scoreProfile(seqs[3], motif,bg)
    expect_equal(c(length(mh[[1]]),length(mh[[2]])),c(10,10))
    expect_equal(mh[[1]],c(-23.0, NaN,NaN,NaN,NaN,NaN,NaN,-14.7,-21.3,-14.8))

    # Check with aaa repeat motif
    name="x8.tab"
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))

    profile=scoreProfile(seqs[2],motif,bg)
    expect_equal(length(profile[[2]]),length(profile[[1]]))
    expect_equal(length(profile[[1]]),length(seqs[[2]])-ncol(motif)+1)
    # in the forward strand one hit
    expect_equal(which(profile[[1]]==max(profile[[1]])),4)

})

test_that("scoreHistogram", {
    # This part implicitly also tests scoreHistogramSingleSeq

    # Set the significance level
    motifcounterOptions(alpha=0.01,gran=0.1)

    # Load sequences
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load background
    bg=readBackground(seqs,1)


    # Load motif
    motiffile=system.file("extdata","x31.tab",package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    # Provoke an error
    expect_error(scoreHistogram(seqfile,motif,bg)) # not a DNAStringSet

    # error because it is no DNAStringSet
    # This was updated. It used to only accept DNAStringSet,
    # but now, DNAStringSet and DNAString is accepted.
    scoreHistogram(seqs[[1]],motif,bg)

    # check motif shorter than bg order
    bg=readBackground(seqs,2)
    expect_error(scoreHistogram(seqs,as.matrix(motif[, 1]),bg))

    # Reload background
    bg=readBackground(seqs,1)

    # Check short sequence must yield only zeros
    dp=scoreDist(as.matrix(motif),bg)
    expect_equal(scoreHistogram(
                Biostrings::DNAStringSet(""),motif,bg)$dist,
                rep(0,length(dp[[1]])))
    expect_equal(scoreHistogram(
                generateDNAStringSet(ncol(motif)-1,bg),motif,bg)$dist,
                rep(0,length(dp[[1]])))

    # scoreHistogram same range as scoreDist
    dp=scoreDist(as.matrix(motif),bg)
    seqs=generateDNAStringSet(10,bg)
    sh=scoreHistogram(seqs,motif,bg)

    expect_equal(sh[[1]],dp[[1]])

    # Check total number of scores
    expect_equal(sum(sh$dist),10-ncol(motif)+1)

    # Insert an 'N' must yield zero observations
    seq=Biostrings::replaceLetterAt(seqs[[1]],3,"N")
    sh=scoreHistogram(Biostrings::DNAStringSet(seq),motif,bg)
    expect_equal(sh$dist,rep(0,length(dp[[1]])))
    expect_equal(sh[[1]],dp[[1]])

    # Test set of sequences
    seqs=generateDNAStringSet(rep(30,100),bg)

    sh=scoreHistogram(seqs,motif,bg)
    expect_equal(sum(sh$dist),(30-ncol(motif)+1)*100) # number of observed

})

test_that("scoreThreshold", {

    # Set the significance level and granularity for the score computation
    motifcounterOptions(alpha=0.01,gran=0.1)

    # Load sequences
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load background
    bg=readBackground(seqs,1)

    # Load  motif
    motiffile=system.file("extdata","x31.tab",package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    # check motif shorter than bg order
    bg=readBackground(seqs,2)
    expect_error(scoreThreshold(as.matrix(motif[, 1]),bg))

    # Reload background
    bg=readBackground(seqs,1)

    # Test plausibility of threshold
    sth=scoreThreshold(motif,bg)
    expect_true(sth$alpha<=0.01)

    # Use different sequence
    seqfile=system.file("extdata","test2.fa", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load background
    bg=readBackground(seqs,0)

    # Load motif
    name="x8.tab"
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))

    # try too stringent
    motifcounterOptions(alpha=0.001,gran=0.1)
    expect_error(scoreThreshold(motif,bg))

    # retry with stringent threshold
    motifcounterOptions(alpha=0.01,gran=0.1)
    scoreThreshold(motif,bg)
})
