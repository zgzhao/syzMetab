## NOTE: some useful functions

.ordnames <- function(x, prefix.fac, n=2) {
    nx <- floor(log10(max(x)))
    nx <- max(c(nx + 1, n))
    rx <- x + 10 ^ nx
    rx <- sub("^1", "", rx)
    paste0(prefix.fac[1], rx)
}

.xlapply <- function(lst, FUN, mc.cores, ...) {
    if(length(lst) < 2) return(lapply(lst, FUN=FUN, ...))
    if(missing(mc.cores)) mc.cores <- detectCores() - 1
    mc.cores <- min(length(lst), mc.cores)
    mclapply(lst, FUN=FUN, mc.cores=mc.cores, ...)
}

.setGeneNames <- function(genes) {
    genes <- unique(unlist(genes))
    ss <- genes %in% c("auto", "HLINK")
    nn <- sum(ss)
    if(nn < 1) return(genes)
    genex <- NULL
    for(i in 1:nn) {
        x <- sample(c(LETTERS, letters, 0:9), 8)
        x <- paste0(x, collapse="")
        genex <- c(genex, x)
    }
    genex <- paste0("auto", genex)
    genes[ss] <- genex
    sort(genes)
}
## CONSTANTS
.nullSTcuts <- new("stcuts", edges=list(), genes=list())
