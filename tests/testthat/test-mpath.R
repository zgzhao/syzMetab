test_that("KDataSet and ReactionList work", {
    pp <- make_kdset("ath00010")
    expect_s4_class(pp, "KDataSet")
    expect_s3_class(pp@reactions, "ReactionList")

    check1 <- sapply(pp@reactions, FUN=function(aa){
        sapply(aa, FUN=function(bb) all(is.vector(bb)))
    })
    expect_equal(all(check1), TRUE)
    check2 <- sapply(pp@reactions, FUN=function(aa){
        all(is.character(aa[["substrate"]]),
            is.character(aa[["product"]]),
            is.character(aa[["gene"]]),
            is.character(aa[["name"]]),
            is.logical(aa[["reversible"]]))
    })
    expect_equal(all(check2), TRUE)
})
