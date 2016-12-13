context("setup")

test_that("Setup check", {
    expect_error(motifcounterOptions(-1))
    expect_error(motifcounterOptions(3))
    expect_warning(motifcounterOptions(.4))
    expect_error(motifcounterOptions(.01,-1))
})
