test_that("class checking works", {
    ll <- list(1:10)
    expect_equal(is.xgraph(ll), FALSE)
    pp <- make_kdset("ath00010")
    gi <- make_empty_graph(10)
    gm <- make_mgraph(pp)
    gg <- make_ggraph(pp)
    gr <- make_rgraph(pp)
    expect_s4_class(pp, "KDataSet")
    expect_s3_class(gi, "igraph")
    expect_s3_class(gm, "mgraph")
    expect_s3_class(gg, "ggraph")
    expect_s3_class(gr, "rgraph")

    gm <- make_mgraph("ath00010")
    gg <- make_ggraph("ath00010")
    gr <- make_rgraph("ath00010")
    expect_s3_class(gm, "mgraph")
    expect_s3_class(gg, "ggraph")
    expect_s3_class(gr, "rgraph")

    chem1 <- c("C00031", "C00221", "C00267", "C01172", "C01451", "C06186")
    chem2 <- "C00022"
    gb <- make_bgraph(gm, s=chem1, t=chem2)
    expect_s3_class(gb, "bgraph")

    xtest <- function(cc, ...) sapply(list(...), function(x) is(x, cc))
    ## S4 methods use "is" fucntion
    ## NOTE: graph list is not "list"
    expect_equal(xtest("list", gi, gm, gg, gr), c(F, F, F, F))
    expect_equal(xtest("igraph", gi, gm, gg, gr), c(T, F, F, F))
    expect_equal(xtest("mgraph", gi, gm, gg, gr), c(F, T, F, F))
    expect_equal(xtest("ggraph", gi, gm, gg, gr), c(F, F, T, F))
    expect_equal(xtest("rgraph", gi, gm, gg, gr), c(F, F, F, T))
    ## NOTE: ensure xgraph methods do not apply to igraph object
    expect_equal(xtest("xgraph", gi, gm, gg, gr), c(F, T, T, T))

    ## S3 methods and manual checking
    xtest <- function(FUN, ...) sapply(list(...), function(x) FUN(x))
    ## NOTE: ensure xgraph object can use igraph S3 methods
    expect_equal(xtest(is.igraph, gi, gm, gg, gr), c(T, T, T, T))
    expect_equal(xtest(is.mgraph, gi, gm, gg, gr), c(F, T, F, F))
    expect_equal(xtest(is.ggraph, gi, gm, gg, gr), c(F, F, T, F))
    expect_equal(xtest(is.rgraph, gi, gm, gg, gr), c(F, F, F, T))
    expect_equal(xtest(is.xgraph, gi, gm, gg, gr), c(F, T, T, T))

    ## ReactionSet
    rr <- make_rset(pp)
    expect_s4_class(rr, "ReactionSet")
    rx <- Reactions(pp)
    rr <- make_rset(rx)
    expect_s4_class(rr, "ReactionSet")
    class(rx) <- "list"
    rr <- make_rset(rx)
    expect_s4_class(rr, "ReactionSet")
    rr <- make_rset(c("ko00010", "ko00020"))
    expect_s4_class(rr, "ReactionSet")
})
