context("Motif tests")

test_that("motif validity",
{
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
    motifValid(normalizeMotif(motif))
    
    # provoke wrong nomalization
    motif[1,1]=.1
    expect_error(motifValid(motif))
    motifValid(normalizeMotif(motif))
})
