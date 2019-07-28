context("Motif tests")

test_that("motif validity",
{
    # Read motif
    motiffile=system.file("extdata","x1.tab", package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))
    motifValid(motif)

    # not a matrix
    expect_error(motifValid("hallo"))

    #wrong dimensionality
    expect_error(motifValid(t(motif)))

    # provoke wrong nomalization
    motif[1,1]=0
    expect_error(motifValid(motif))
    
    # renormalize to solve the issue
    motifValid(normalizeMotif(motif))
    
    # provoke wrong nomalization
    motif[1,1]=.1
    expect_error(motifValid(motif))

    # renormalize to solve the issue
    motifValid(normalizeMotif(motif))

    # test reverse complement
    expect_equal(as.vector(motif), as.vector(revcompMotif(revcompMotif(motif))))
})


test_that("position weight matrix",
{
    seqfile = system.file("extdata", "seq.fasta", package="motifcounter")
    seqs=Biostrings::readDNAStringSet(seqfile)

    motiffile=system.file("extdata", "x1.tab", package="motifcounter")
    motif=t(as.matrix(read.table(motiffile)))

    bg = readBackground(seqs, 0)
    pwm = getPositionWeights(motif, bg)

    expect_equal(dim(pwm), dim(motif))

    expect_equal(pwm, matrix(c(-55, -55, -55, 14,
                               -55, -55,  14, -55,
                               -55, -55,  14, -55,
                               -55, -55, -55,  14,
                               -55,  14, -55, -55,
                               -55,  14, -55, -55), ncol=6, nrow=4))

    bg = readBackground(seqs, 1)
    pwm = getPositionWeights(motif, bg)

    expect_equal(dim(pwm), c(16, ncol(motif) -1))

    expect_equal(pwm[1,1], -111)
    expect_equal(pwm[15,1], 28)
    expect_equal(pwm[16,1], -42)
    expect_equal(pwm[16,5], -56)
})
