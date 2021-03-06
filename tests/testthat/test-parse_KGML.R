test_that("test KGML parsing", {
    pp1 <- make_kdset("ko00010")
    ## entries
    entx <- pp1@entries[["52"]]
    expect_identical(sort(entx[["name"]]), c("K00161", "K00162", "K00163"))
    expect_identical(entx[["type"]], "ortholog")
    expect_identical(entx[["reaction"]], "R00014")
    entx <- pp1@reactions[["52"]]
    ## reactions
    expect_identical(entx[["name"]], "R00014")
    expect_identical(entx[["reversible"]], FALSE)
    expect_identical(sort(entx[["substrate"]]), c("C00022", "C00068"))
    expect_identical(entx[["product"]], "C05125")
    expect_identical(sort(entx[["gene"]]), c("K00161", "K00162", "K00163"))
    ## fit organism: more reactions from KO-path
    pp2 <- make_kdset("ath00010")
    pp3 <- setorg_path(pp1, "ath")
    expect_identical(length(pp2@reactions) <= length(pp3@reactions), TRUE)
})
