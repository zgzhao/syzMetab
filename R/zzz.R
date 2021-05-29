## library(KNetStats)
## kegdata <- keggs
## kegdata$parents <- apply(kegdata[, 4:6], 1, FUN = function(x) setdiff(x, NA))
## kegdata <- kegdata[, -(4:6)]
## kegdata[1, ]
## save(kegdata, file="../data/kegdata.rda")
