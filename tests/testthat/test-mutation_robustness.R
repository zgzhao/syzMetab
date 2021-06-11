test_that("mutational robustness", {
    s <- c("C00221", "C00267", "C00103")
    p <- "C00022"
    gg <- make_mgraph("ath00010")
    aa <- make_bgraph(gg, s, p)
    bb <- all_stcuts(gg, s, p)
    rbs1 <- robust(aa)
    rbs2 <- robust(bb)
    expect_equal(rbs1, rbs2)
})
