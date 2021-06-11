#' Calculate functional probability of a gene cutset
#'
#'  Note that a "functional" cutset of genes means all genes in it are mutated. If containing overlap genes, the cut sets must be grouped and joint probability is calculated.
#' @title functional probability of s-t gene cutsets
#' @param g.sets list, each item is a named gene vector and the names (set names) are used for grouping.
#' @param s.groups list, each item if a vector containg set names.
#' @param GMR numeric, gene mutation rate (default 0.01).
#' @return numeric vector, functional probability of each gene group.
#' @author ZG Zhao
#' @export
prob_stcuts <- function(g.sets, s.groups, GMR=0.01) {
    ppx <- lapply(s.groups, FUN=function(ffs){
        glist <- g.sets[ffs]
        nfs <- length(glist)
        genes1 <- glist[[1]]
        pxx <- calGroupProb(genes1, GMR=GMR)
        if(nfs < 2) return(pxx)
        for(i in 2:nfs) {
            genes2 <- genes1
            genes1 <- glist[[i]]
            pxx <- calGroupProb(genes1, genes2, pxx, GMR)
        }
        pxx
    })
    ppx <- unlist(ppx)
    names(ppx) <- sapply(s.groups, paste0, collapse="-")
    ppx
}

calGroupProb <- function(genes1, genes2=NULL, p=1, GMR=0.01) {
    nn <- length(genes1)
    if(p == 1) return(GMR ^ nn)
    genex <- union(genes1, genes2)
    n1 <- length(genes2)
    n2 <- length(genex)
    p + GMR ^ nn - GMR ^ n2
}
