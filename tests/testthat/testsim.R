context("Simulation")

test_that("simulate DNA", {
    deleteBackgroundForSampling()
    expect_error(generateDNAString(10)) # load background first

    seqfile=system.file("extdata","seq.fasta", package="motifcounter")
    readBackgroundForSampling(seqfile,1)
    generateDNAString(100) # after background loading everything should work

    expect_error(generateDNAString(0)) # positive length

    readBackgroundForSampling(seqfile,2)
    expect_error(generateDNAString(1)) # len must be at least 2
})
