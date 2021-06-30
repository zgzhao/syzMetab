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
        pxx <- NULL
        for(i in 1:length(glist)) {
            genes <- glist[[i]]
            pxx <- calGroupProb(genes, pxx, GMR)
        }
        pxx[["pval"]]
    })
    ppx <- unlist(ppx)
    names(ppx) <- sapply(s.groups, paste0, collapse="-")
    ppx
}

calGroupProb <- function(genes, plast=NULL, GMR=0.01) {
    genes <- unique(genes)
    pp <- GMR ^ length(genes)
    if(is.empty(plast)) return(list(pval=pp, genes=genes))
    genex <- union(genes, plast$genes)
    pp <- plast$pval + pp - GMR ^ length(genex)
    list(pval=pp, genes=genex)
}
