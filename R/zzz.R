## NOTE: some useful functions

.ordnames <- function(x, prefix.fac, n=3) {
    nx <- ceiling(log10(max(x)))
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

## return a graph without nodes and edges
.emptyGraph <- function(g){
    delete.vertices(g, vnames(g))
}
