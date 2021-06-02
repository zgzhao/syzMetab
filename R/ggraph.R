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
    for(ndx in names(rlist)) {
        rx <- rlist[[ndx]]
        gns <- rx[["gene"]]
        pp <- rx[["product"]]
        gnp <- lapply(rlist, FUN=function(rr){
            ss <- rr[["substrate"]]
            if(any(ss %in% pp)) return(rr[["gene"]])
            else return(NULL)
        })
        gnp <- unlist(gnp)
        gnp <- setdiff(gnp, gns)
        g <- xaddEdges(g, gns, gnp)
    }
    attr(g, "reactions") <- object
    class(g) <- c("ggraph", class(g))
    g
})
setMethod("make_ggraph", "character", function(object, d.path="KEGG"){
    rsobj <- rset_from_kos(object, d.path)
    make_ggraph(rsobj)
})
setMethod("make_ggraph", "keggPATH", function(object){
    rsobj <- as_rset(object)
    make_ggraph(rsobj)
})

setMethod("make_ggraph", "ReactionList", function(object, org){
    rsobj <- as_rset(object, org)
    make_ggraph(rsobj)
})
