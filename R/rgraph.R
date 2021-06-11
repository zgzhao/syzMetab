#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
# File: rgraph.R
# Description: reaction graph
# AUTHOR: ZG Zhao; zgzhao@foxmail.com
# 2021-05-31 17:36:10

## refer to make_mgraph for manul

#' @export
setGeneric("make_rgraph", function(object, ...) standardGeneric("make_rgraph"))
setMethod("make_rgraph", "ReactionSet", function(object){
    vv <- rnames(object)
    g <- make_empty_graph(n=length(vv))
    V(g)$name <- vv
    rlist <- Reactions(object)
    for(aa in vv) {
        rx <- rlist[[aa]]
        pp <- rx[["product"]]
        bb <- sapply(vv, FUN=function(x){
            if(x == aa) return(FALSE)
            ss <- rlist[[x]][["substrate"]]
            any(ss %in% pp)
        })
        g <- add.edges(g, expand.grid(aa, vv[bb]))
    }
    attr(g, "reactions") <- object
    class(g) <- c("rgraph", class(g))
    g
})

setMethod("make_rgraph", "character", function(object, d.path="KEGG"){
    rsobj <- make_rset(object, d.path)
    make_rgraph(rsobj)
})
setMethod("make_rgraph", "keggPATH", function(object){
    rsobj <- make_rset(object)
    make_rgraph(rsobj)
})

setMethod("make_rgraph", "ReactionList", function(object, org){
    rsobj <- make_rset(object, org)
    make_rgraph(rsobj)
})
