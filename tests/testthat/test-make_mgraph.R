test_that("make mgraph from various objects", {
    gg1 <- make_mgraph("ath00010")
    kk <- kmeta_from_ko("ath00010")
    gg2 <- make_mgraph(kk)
    kr <- kk@reactions
    gg3 <- make_mgraph(kr, "ath")
    expect_equal(isomorphic(gg1, gg2), TRUE)
    expect_equal(isomorphic(gg1, gg3), TRUE)
})
