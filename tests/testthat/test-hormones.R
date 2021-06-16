test_that("robust in practice", {
    ## kos <- read.csv("KEGG/KO-hormones.csv", row.names = 1)
    ## chems <- read.csv("KEGG/CHEM-hormones.csv")
    ## hormones <- rownames(kos)
    ## for (hh in hormones) {
    ##     kox <- strsplit(kos[hh, "ko"], ";")[[1]]
    ##     chem1 <- chems[chems$HORMONE == hh & chems$TYPE == "S", "CPD"]
    ##     chem2 <- chems[chems$HORMONE == hh & chems$TYPE == "E", "CPD"]
    ##     if(length(chem1) < 1) next
    ##     gg <- make_mgraph(kox)
    ##     ffs <- file.path("KEGG/XLNK", paste0(kox, ".csv"))
    ##     for(ff in ffs) {
    ##         if(file.exists(ff)) {
    ##             rr <- read.csv(ff)
    ##             gg <- add.reactions(gg, rr)
    ##         }
    ##     }
    ##     gx <- setorg_graph(gg, "ath")
    ##     gb <- make_bgraph(gx, chem1, chem2)
    ##     cat(hh, robust(gb), "\n")
    ## }
})
