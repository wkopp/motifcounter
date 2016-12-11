context("Score distribution")

test_that("scorethreshold", {
    # Set the significance level and granularity for the score computation
    motifcounterOption(alpha=0.01,gran=0.1)
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    # Load the order-1 background model from the DNA sequence
    bg=readBackground(seqs,1)

    motiffile=system.file("extdata","x31.tab",package="motifcounter")

    # Load the motif from the motiffile
    motif=t(as.matrix(read.table(motiffile)))

    sth=scoreThreshold(motif,bg)
    expect_true(sth$alpha<=0.01)

    name="x8.tab"
    seqfile=system.file("extdata","test2.fa", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    file=system.file("extdata",name, package="motifcounter")
    motif=t(as.matrix(read.table(file)))

    bg=readBackground(seqs,0)

    motifcounterOption(alpha=0.001,gran=0.1)
    expect_error(scoreThreshold(motif,bg)) # too stringent
    motifcounterOption(alpha=0.01,gran=0.1)
    scoreThreshold(motif,bg)
})

test_that("scorehistogram", {
    # Set the significance level and granularity for the score computation
    motifcounterOption(alpha=0.01,gran=0.1)

    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    # Load the order-1 background model from the DNA sequence
    bg=readBackground(seqs,1)


    seq=generateDNAString(10,bg)

    motiffile=system.file("extdata","x31.tab",package="motifcounter")

    # Load the motif from the motiffile
    motif=t(as.matrix(read.table(motiffile)))

    expect_error(scoreHistogram(seqfile,motif,bg)) # not a DNAStringSet

    # scoreHistogram same range as scoreDist
    dp=scoreDist(as.matrix(motif),bg)
    sh=scoreHistogram(Biostrings::DNAStringSet(seq),motif,bg)

    expect_equal(sh[[1]],dp[[1]])

    # expect non-zero entries
    expect_equal(sum(sh$frequency),10-ncol(motif)+1)

    seq=Biostrings::replaceLetterAt(seq,3,"N")
    sh=scoreHistogram(Biostrings::DNAStringSet(seq),motif,bg)
    expect_equal(sum(sh$frequency),0)

    seqs=generateDNAStringSet(rep(30,100),bg)

    sh=scoreHistogram(seqs,motif,bg)
    expect_equal(sum(sh$frequency),(30-ncol(motif)+1)*100) # number of observed

    bg=readBackground(seqs,2)
    expect_error(scoreHistogram(seqs[[1]][1],motif,bg))
})

test_that("scoresequence", {
    # Set the significance level and granularity for the score computation
    motifcounterOption(alpha=0.01,gran=0.1)

    # Load the order-1 background model from the DNA sequence
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    expect_error(scoreSequence(seqs,motif,bg))
    bg=readBackground(seqs,1)


    seq=generateDNAString(10,bg)


    motiffile=system.file("extdata","x31.tab",package="motifcounter")

    # Load the motif from the motiffile
    motif=t(as.matrix(read.table(motiffile)))

    # Compute the score distribution
    scores=scoreSequence(seq,motif,bg)
    # expect error with too short sequence
    expect_error(scoreSequence(generateDNAString(ncol(motif)-1,bg),motif,bg))

    # same number of score on both strands because its a palindrome
    expect_equal(length(scores$fscores),length(scores$rscores))

    # Nscores= Nseq-Nmotif+1
    expect_equal(length(scores$fscores),10-ncol(motif)+1)

    seq=Biostrings::replaceLetterAt(seq,3,"N")
    scores=scoreSequence(seq,motif,bg)
    expect_equal(scores$fscores,scores$rscores) # both should be equal 
})

test_that("score correctness", {
    alpha=0.01
    motifcounterOption(alpha)
    motifname="x1.tab"
    motiffile=system.file("extdata",motifname, package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    # 1. test: with motif of length one

    bg=readBackground(seqs, 0)

    # compute the score distribution for this case in R
    s=round((log(motif[,1])-log(bg$station))*10)
    srange=seq(min(s),max(s))
    p=rep(0,length(srange))
    for (i in 1:length(srange)){
        p[i]=sum(bg$station[which(s==srange[i])])
    }

    dp=scoreDist(as.matrix(motif[,1]),bg)

    expect_equal(sum(dp[[2]]),1)
    expect_equal(length(srange),length(dp[[1]]))
    expect_equal(srange,round(dp[[1]]*10))
    expect_equal(p,dp[[2]])


    # test whether the range is equally long
    # test whether the zero entries in the score
    # distribution overlap perfectly
    # test with the stationary probabilities of the background model only

    for (m in seq(0,3)) {
        bg=readBackground(seqs, m)

        if (m==0) {
            #simluate score distribution
            sims=scoreDistEmpirical(as.matrix(motif[,1]), bg,1,1600)
            #compute exact score distribution
            dp=scoreDist(as.matrix(motif[,1]),bg)
            bf=scoreDistBf(as.matrix(motif[,1]),bg)

        } else {
            #simluate score distribution
            sims=scoreDistEmpirical(as.matrix(motif[,1:m]),bg,m,1600)
            #compute exact score distribution
            dp=scoreDist(as.matrix(motif[,1:m]),bg)
            bf=scoreDistBf(as.matrix(motif[,1:m]),bg)
        }

        expect_equal(sum(dp[[2]]),1)
        expect_equal(length(bf[[1]]),length(dp[[1]]))
        expect_equal(bf[[1]],dp[[1]])
        expect_equal(length(bf[[2]]),length(dp[[2]]))
        expect_equal(bf[[2]],dp[[2]])
        expect_equal(length(sims[[1]]),length(dp[[1]]))
        expect_equal(sims[[1]],dp[[1]])
        expect_true(all(!xor(sims[[2]]==0,dp[[2]]==0)))
    }

    # test whether the range is equally long
    # test whether the zero entries in the score distribution overlap perfectly
    # test with the stationary and the transition probabilities
    for (m in seq(1,3)) {
        bg=readBackground(seqs, m)


        #simluate score distribution
        #compute exact score distribution
        dp=scoreDist(motif[,1:(m+1)],bg)
        bf=scoreDistBf(motif[,1:(m+1)],bg)


        expect_equal(sum(dp[[2]]),1)
        expect_equal(length(bf[[1]]),length(dp[[1]]))
        expect_equal(bf[[1]],dp[[1]])
        expect_equal(length(bf[[2]]),length(dp[[2]]))
        expect_equal(bf[[2]],dp[[2]])

        sims=scoreDistEmpirical(motif[,1:(m+1)],bg,m+1,3000)
        expect_equal(length(sims[[1]]),length(dp[[1]]))
        expect_equal(sims[[1]],dp[[1]])
        expect_true(all(!xor(sims[[2]]==0,dp[[2]]==0)))

    }
})

test_that("score correctness with jaspar motifs", {
    library(MotifDb)
    motifs=as.list(query(query(MotifDb,"hsapiens"),"JASPAR_CORE"))

    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    bg=readBackground(seqs, 1)

    for ( motif in motifs ) {
        motif=normalizeMotif(motif)
        if (ncol(motif)>4) {
            motif=motif[, 1:4]
        }
        dp=scoreDist(as.matrix(motif),bg)
        bf=scoreDistBf(as.matrix(motif),bg)

        expect_equal(dp$probability, bf$probability)
        expect_equal(sum(dp$probability),1)

        expect_true(dp$probability[1]>0)
        expect_true(dp$probability[length(dp[[2]])]>0)
    }

    logp1=round(log(motif[,1]/bg$station)*10)
    logp2=round(log(rep(motif[,2],4)/bg$trans)*10)
    score=as.vector(matrix(logp2,4,4)+t(matrix(rep(logp1,4),4,4)))
    p1=bg$station
    p12=as.vector(matrix(rep(bg$station,4)*bg$trans,4,4))
})
