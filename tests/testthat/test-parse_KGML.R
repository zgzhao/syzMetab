test_that("test KGML parsing", {
    kp1 <- kmeta_from_ko("ko00010")
    ## entries
    entx <- kp1@entries[["52"]]
    expect_identical(sort(entx[["name"]]), c("K00161", "K00162", "K00163"))
    expect_identical(entx[["type"]], "ortholog")
    expect_identical(entx[["reaction"]], "R00014")
    entx <- kp1@reactions[["52"]]
    ## reactions
    expect_identical(entx[["name"]], "R00014")
    expect_identical(entx[["reversible"]], FALSE)
    expect_identical(sort(entx[["substrate"]]), c("C00022", "C00068"))
    expect_identical(entx[["product"]], "C05125")
    expect_identical(sort(entx[["gene"]]), c("K00161", "K00162", "K00163"))
    ## fit organism
    kp2 <- kmeta_from_ko("hsa00010")
    kp3 <- kmeta_x_org(kp1, "hsa")
    expect_identical(kp2@reactions, kp3@reactions)
})
