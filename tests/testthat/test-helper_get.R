test_that("test getXXX helpers", {
    pp <- make_kdset("ath00010")
    gg <- make_mgraph(pp)
    expect_equal(Organism(pp), "ath")
    expect_equal(Organism(gg), "ath")
    expect_equal(Chemicals(pp), Chemicals(gg))
    expect_equal(Genes(pp), Genes(gg))
})
