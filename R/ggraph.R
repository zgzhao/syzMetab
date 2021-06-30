#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
# File: ggraph.R
# Description: handling gene network
# AUTHOR: ZG Zhao; zgzhao@foxmail.com
# 2021-05-31 16:16:37

## refer to make_mgraph for manul

#' @export
setGeneric("make_ggraph", function(object, ...) standardGeneric("make_ggraph"))
setMethod("make_ggraph", "ReactionSet", function(object){
    vv <- Genes(object)
    g <- make_empty_graph(n=length(vv))
    V(g)$name <- vv

    rlist <- Reactions(object)
    for(rx in rlist) {
        pp <- rx[["product"]]
        gene1 <- rx[["gene"]]
        gene2 <- lapply(rlist, FUN=function(rr){
            ss <- rr[["substrate"]]
            if(any(ss %in% pp)) return(rr[["gene"]])
            else return(NULL)
        })
        gene2 <- unlist(gene2)
        gene2 <- setdiff(gene2, gene1)
        if(is.empty(gene2)) next
        dd <- expand.grid(gene1, gene2)
        g <- add.edges(g, dd)
    }
    attr(g, "reactions") <- object
    class(g) <- c("ggraph", class(g))
    g
})
setMethod("make_ggraph", "character", function(object, d.path="KEGG"){
    rsobj <- make_rset(object, d.path)
    make_ggraph(rsobj)
})
setMethod("make_ggraph", "KDataSet", function(object){
    rsobj <- make_rset(object)
    make_ggraph(rsobj)
})

setMethod("make_ggraph", "ReactionList", function(object, org){
    rsobj <- make_rset(object, org)
    make_ggraph(rsobj)
})
