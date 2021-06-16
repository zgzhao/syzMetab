test_that("mutational robustness", {
    s <- c("C00221", "C00267", "C00103")
    t <- "C00022"

    ## setorg and stcuts
    gg <- make_mgraph("ko00010")
    aa <- setorg_graph(gg, "ath")
    bb <- make_stcuts(gg, s, t)
    rbs1 <- robust(make_stcuts(aa, s, t))
    rbs2 <- robust(setorg_stcuts(bb, "ath"))
    expect_equal(rbs1, rbs2)

    ## bgraph and stcuts
    gg <- make_mgraph("ath00010")
    rbs3 <- robust(make_bgraph(gg, s, t))
    rbs4 <- robust(make_stcuts(gg, s, t))
    expect_equal(rbs2, rbs3)
    expect_equal(rbs3, rbs4)

    ## with or without extracting simple paths
    gx1 <- subset_graph(gg, s, t, clean=FALSE)
    gx2 <- subset_graph(gg, s, t, clean=TRUE)
    expect_identical(ecount(gx1) > ecount(gx2), TRUE)
    gb1 <- make_bgraph(gx1, s, t)
    gb2 <- make_bgraph(gx2, s, t)
    rbs5 <- robust(gb1)
    rbs6 <- robust(gb2)

    expect_equal(rbs4, rbs5)
    expect_equal(rbs5, rbs6)

})
