test_that("make mgraph from various objects", {
    pp <- make_mpath("ath00010")
    gg1 <- make_mgraph(pp)
    gg2 <- make_mgraph("ath00010")
    rr <- pp@reactions
    gg3 <- make_mgraph(rr, "ath")
    expect_equal(isomorphic(gg1, gg2), TRUE)
    expect_equal(isomorphic(gg1, gg3), TRUE)
})
