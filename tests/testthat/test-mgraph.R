test_that("class checking works", {
    ## KDataSet
    pp1 <- make_kdset("ko00010")
    pp2 <- make_kdset("ko00020")
    ## ReactionList/list
    rr1 <- Reactions(pp1)
    ## ReactionSet
    rr2 <- Reactions(pp1)
    class(rr2) <- "list"
    ## character (ko)
    gm1 <- make_mgraph(pp1)
    gm2 <- make_mgraph(rr1)
    gm3 <- make_mgraph(rr2)
    gm4 <- make_mgraph("ko00010")
    expect_equal(isomorphic(gm1, gm2), TRUE)
    expect_equal(isomorphic(gm2, gm3), TRUE)
    expect_equal(isomorphic(gm3, gm4), TRUE)

    ## list
    rr3 <- c(Reactions(pp1), Reactions(pp2))
    class(rr3) <- "list"
    gm5 <- make_mgraph(rr3)
    gm6 <- make_mgraph(c("ko00010", "ko00020"))
    expect_equal(isomorphic(gm5, gm6), TRUE)

    ## check reaction lists
    rrx <- Reactions(gm6)
    check1 <- sapply(rrx, FUN=function(aa){
        all(is.vector(aa[["substrate"]]),
            is.vector(aa[["product"]]),
            is.vector(aa[["gene"]]))
    })
    expect_equal(all(check1), TRUE)
    check2 <- sapply(rrx, FUN=function(aa){
        all(is.character(aa[["substrate"]]),
            is.character(aa[["product"]]),
            is.character(aa[["gene"]]))
    })
    expect_equal(all(check2), TRUE)
    check3 <- sapply(rrx, FUN=function(aa){
        xx <- aa[c("substrate", "product")]
        length(unlist(xx)) == 2
    })
    expect_equal(all(check3), TRUE)

})
