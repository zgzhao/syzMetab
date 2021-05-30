test_that("test getXXX helpers", {
    kk <- kmeta_from_ko("ath00010")
    gg <- make_mgraph(kk)
    expect_equal(Organism(kk), "ath")
    expect_equal(Organism(gg), "ath")
    expect_equal(Chemicals(kk), Chemicals(gg))
    expect_equal(Genes(kk), Genes(gg))
    ## NOTE: reactions in KEGGmeta and mgraph are different!
    ## expect_equal(Reactions(kk), Reactions(gg))
})
