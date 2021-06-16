test_that("test gene name settings", {
    gm1 <- make_mgraph("ko00010")
    gm2 <- add.reactions(gm1, "C00186", "C00074", "auto")
    expect_equal(any(grepl("^auto", Genes(gm2))), TRUE)
    expect_equal(any(grepl("^auto$", Genes(gm2))), FALSE)

    gm3 <- setorg_graph(gm2, "ath")
    expect_equal(any(grepl("^absent$", Genes(gm3))), TRUE)
    expect_equal(any(grepl("^auto", Genes(gm3))), TRUE)
    expect_equal(any(grepl("^auto$", Genes(gm3))), FALSE)

    s <- c("C00221", "C00267", "C00103")
    t <- "C00022"
    gx <- subset_graph(gm3, s, t, clean = FALSE)
    expect_equal(any(grepl("^absent$", Genes(gx))), FALSE)
    expect_equal(any(grepl("^auto", Genes(gx))), TRUE)
    expect_equal(any(grepl("^auto$", Genes(gx))), FALSE)
})
