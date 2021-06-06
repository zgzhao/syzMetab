test_that("mutational robustness", {
    s <- c("C00221", "C00267", "C00103")
    p <- "C00022"
    gg <- make_mgraph("ath00010")
    ## set chemicals (clean network) should not alter basic properpties of the network
    gx1 <- xgraph_setchems(gg, s, p, clean=FALSE)
    gx2 <- xgraph_setchems(gg, s, p, clean=TRUE)
    ## gx1 should be larger than gx2
    expect_equal(vcount(gx1) > vcount(gx2), TRUE)
    expect_equal(ecount(gx1) > ecount(gx2), TRUE)
    ##
    rbs1 <- robust(gx1, s, p, GMR=0.001)
    rbs2 <- robust(gx2, s, p, GMR=0.001)
    expect_equal(rbs1, rbs2)
})
