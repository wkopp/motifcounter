context("setup")

test_that("Setup check", {
    expect_error(motifcounterOption(-1))
    expect_error(motifcounterOption(3))
    expect_warning(motifcounterOption(.4))
    expect_error(motifcounterOption(.01,-1))
})
