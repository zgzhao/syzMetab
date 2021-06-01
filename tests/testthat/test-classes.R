test_that("class checking works", {
    ll <- list(1:10)
    expect_equal(is.xgraph(ll), FALSE)

    pp <- make_mpath("ath00010")
    gi <- make_empty_graph(10)
    gm <- make_mgraph(pp)
    gg <- make_ggraph(pp)
    gr <- make_rgraph(pp)

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
})
