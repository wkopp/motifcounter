context("setup")

test_that("motifcounterOptions", {
    # negative alpha
    expect_error(motifcounterOptions(-1))
    # too large alpha
    expect_error(motifcounterOptions(3))
    # warning with large alpha
    expect_warning(motifcounterOptions(.4))
    # negative granularity
    expect_error(motifcounterOptions(.01,-1))
    # zero granularity
    expect_error(motifcounterOptions(.01,0.0))

    # retrieve alpha
    alpha=0.01
    motifcounterOptions(alpha)
    expect_equal(alpha,sigLevel())
})
