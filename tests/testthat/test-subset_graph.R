test_that("subset_graph", {
    chems <- read.csv("KEGG/CHEM-hormones.csv")
    kos <- read.csv("KEGG/KO-hormones.csv", row.names = 1)
    hh <- "GA"
    kox <- strsplit(kos[hh, "ko"], ";")[[1]]
    chem1 <- chems[chems$HORMONE == hh & chems$TYPE == "S", "CPD"]
    chem2 <- chems[chems$HORMONE == hh & chems$TYPE == "E", "CPD"]
    gg <- make_mgraph(kox)
    ffs <- file.path("KEGG/XLNK", paste0(kox, ".csv"))
    for(ff in ffs) {
        if(file.exists(ff)) {
            rr <- read.csv(ff)
            gg <- add.reactions(gg, rr)
        }
    }
    gx1 <- setorg_graph(gg, "aly")
    gx2 <- subset_graph(gx1, chem1, chem2, clean=TRUE)
    gx3 <- make_bgraph(gx1, chem1, chem2)
    expect_equal(is.empty(gx2), FALSE)
    expect_equal(is.empty(gx3), FALSE)
})
