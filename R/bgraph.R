#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
# File: bgraph.R
# Description: bipartie graph contains both metabolites and gene sets
# AUTHOR: ZG Zhao; zgzhao@foxmail.com
# 2021-06-10 07:40:25


#' @export
setGeneric("make_bgraph", function(object, s, t, prefix.fac="X", ...) standardGeneric("make_bgraph"))
setMethod("make_bgraph", "NULL", function(object, s, t, prefix.fac) return(NULL))
setMethod("make_bgraph", "stcuts", function(object, prefix.fac){
    gsets <- object@genes
    gnn <- sapply(gsets, length)
    genes <- unique(unlist(gsets))
    n1 <- length(gsets)
    n2 <- length(genes)
    isfac <- c(rep(TRUE, n1), rep(FALSE, n2))
    s.names <- .ordnames(1:n1, prefix.fac)
    names(gsets) <- s.names

    g <- make_bipartite_graph(types=isfac, edges=NULL)
    V(g)$name <- c(s.names, genes)
    for(sx in s.names) g <- add.edges(g, expand.grid(sx, gsets[[sx]]))
    V(g)$isfac <- isfac
    object@genes <- gsets
    attr(g, "stcuts") <- object
    class(g) <- c("bgraph", class(g))
    g
})
setMethod("make_bgraph", "mgraph", function(object, s, t, prefix.fac){
    xcuts <- make_stcuts(object, s, t)
    make_bgraph(xcuts, prefix.fac=prefix.fac)
})
