
context("Markov model")

test_that("markovModel", {
    # inital settings

    alpha=0.01
    gran=0.1
    motifcounterOptions(alpha, gran)

    # Obtain background
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)

    # Obtain motif
    pwmname="x3.tab"
    motiffile=system.file("extdata",pwmname, package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    # overlapping hit probs
    op=probOverlapHit(motif,bg, singlestranded = TRUE)
    #print(getAlpha(op))

    # Dist must sum to one
    # single seq
    val = rep(0, ncol(motif))
    val[1] = 1
    expect_equal(markovModel(op, 0)$dist, val)

    expect_equal(sum(markovModel(op, 1)$dist), 1)
    expect_equal(sum(markovModel(op, 10)$dist), 1)
    expect_equal(sum(markovModel(op, 100)$dist), 1)


    # overlapping hit probs
    op=probOverlapHit(motif,bg, singlestranded = FALSE)

    # Dist must sum to one
    # single seq
    val = rep(0, ncol(motif)*2 + 2)
    val[1] = 1
    expect_equal(markovModel(op, 0)$dist, val)

    expect_equal(sum(markovModel(op, 1)$dist), 1)
    expect_equal(sum(markovModel(op, 10)$dist), 1)
    expect_equal(sum(markovModel(op, 100)$dist), 1)

})

test_that("clumpStartProbSingleStranded", {
    # inital settings

    alpha=0.01
    gran=0.1
    motifcounterOptions(alpha, gran)

    # Obtain background
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)

    # Obtain motif
    pwmname="x32.tab" # palindrom
    motiffile=system.file("extdata",pwmname, package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    # overlapping hit probs
    op=probOverlapHit(motif,bg, singlestranded = TRUE)
    expect_true(computeClumpStartProb(op)<getAlpha(op))

    pwmname="x21.tab" # A-repeat
    motiffile=system.file("extdata",pwmname, package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    # overlapping hit probs
    op=probOverlapHit(motif,bg, singlestranded = TRUE)
    expect_true(computeClumpStartProb(op)<getAlpha(op))

    pwmname="x1.tab" # Non-self-overlapping
    motiffile=system.file("extdata",pwmname, package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    # overlapping hit probs
    op=probOverlapHit(motif,bg, singlestranded = TRUE)
    expect_true(computeClumpStartProb(op)>getAlpha(op))

    pwmnames = c("x32.tab", "x21.tab", "x1.tab")

    for (pwmname in pwmnames) {
        motiffile=system.file("extdata",pwmname, package="motifcounter")
        motif=t(as.matrix(read.table(motiffile)))

        # overlapping hit probs
        op=probOverlapHit(motif,bg, singlestranded = TRUE)
        alpha = getAlpha(op)
        #print(computeClumpStartProb(op))
        op@alpha = computeClumpStartProb(op)
        expect_equal(markovModel(op, 100)$dist[2], alpha)
    }


})


test_that("clumpStartProbDoubleStranded", {
    # inital settings

    alpha=0.01
    gran=0.1
    motifcounterOptions(alpha, gran)

    # Obtain background
    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)
    bg=readBackground(seqs,1)

    # Obtain motif
    pwmname="x32.tab" # palindrom
    motiffile=system.file("extdata",pwmname, package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    # overlapping hit probs
    op=probOverlapHit(motif,bg, singlestranded = FALSE)
    expect_true(computeClumpStartProb(op)<getAlpha(op))

    pwmname="x21.tab" # A-repeat
    motiffile=system.file("extdata",pwmname, package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    # overlapping hit probs
    op=probOverlapHit(motif,bg, singlestranded = FALSE)
    expect_true(computeClumpStartProb(op)<getAlpha(op))

    pwmnames = c("x32.tab", "x21.tab", "x1.tab")

    for (pwmname in pwmnames) {
        motiffile=system.file("extdata",pwmname, package="motifcounter")
        motif=t(as.matrix(read.table(motiffile)))

        # overlapping hit probs
        op=probOverlapHit(motif,bg, singlestranded = FALSE)
        alpha = getAlpha(op)
        #print(computeClumpStartProb(op))
        op@alpha = computeClumpStartProb(op)
        dist = markovModel(op, 100)$dist
        expect_equal(dist[2] + dist[3], 2*alpha)
    }
})
